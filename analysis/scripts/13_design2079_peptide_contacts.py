#!/usr/bin/env python3
"""
analysis/scripts/13_design2079_peptide_contacts.py

Per-residue peptide-contact analysis for the CYCLE 02 hero (design_2079_seq00).
Mirror of 12_design3010_peptide_contacts.py, adapted to the cycle-02 FULL 4-chain
layout. Closes the open thread: did the de novo cycle-02 hero actually READ the
peptide, or is it (like the cycle-03 BAKER-scaffold binders) a peptide-blind
MHC-framework binder?

Cycle-02 layout (full 4-chain target):
    A = HLA-A*02:01 heavy chain (~275 aa)
    B = beta-2 microglobulin   (~100 aa)
    C = peptide RMFPNAPYL      (9 aa)
    D = binder (designed)      (~99 aa)

Metric definitions (match the interface-8A ruler used for the controls and for
c02 metrics.jsonl -- HADRUP_JENKINS_2025 Fig 1B style):
    iface_* = mean over binder<->target residue pairs whose min heavy-atom
              distance <= 8 A, of min(PAE_ij, PAE_ji).
              Decomposed into pep (peptide only), mhc (HLA+b2m), tot (all target).
    iface_pep = inf  <=>  no binder atom within 8 A of any peptide atom.
    ipLDDT    = mean pLDDT over residues participating in any <=8 A interface pair.

Run on pod:
    export TMPDIR=/workspace/.cache/tmp; unset XLA_FLAGS; unset LD_LIBRARY_PATH
    cd /workspace/pipeline
    uv run python analysis/scripts/13_design2079_peptide_contacts.py

Optional overrides:
    --pred-dir results/cycle_02/stage2/af2_designs/predictions
    --design   design_2079_seq00
    --binder-chain D  --mhc-chains A,B
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

PEPTIDE_SEQ = "RMFPNAPYL"

# RMFPNAPYL exposure annotation on HLA-A*02:01, reused verbatim from the 3010 run
PEP_ROLE = {
    1: "N-term",
    2: "ANCHOR buried",
    3: "semi-exposed",
    4: "EXPOSED",
    5: "EXPOSED-specificity",
    6: "EXPOSED",
    7: "EXPOSED",
    8: "EXPOSED-specificity",
    9: "ANCHOR buried",
}
EXPOSED_POS = {4, 5, 6, 7, 8}
ANCHOR_POS = {2, 9}

TIGHT_A = 5.0
CONTACT_A = 8.0

DEFAULT_PRED_DIR = "results/cycle_02/stage2/af2_designs/predictions"
DEFAULT_DESIGN = "design_2079_seq00"

# Reference bands for inline comparison (interface-8A, from phase_D_findings.md)
POS_TOT = (1.27, 1.71)
POS_PEP = (1.46, 2.60)
REF_3010 = dict(tot=2.98, pep=2.10, mhc=3.35)
REF_3054 = dict(tot=1.61, pep=float("inf"))


def find_files(pred_dir: Path, design: str):
    pdbs = sorted(pred_dir.glob(f"*{design}*rank_001*.pdb"))
    if not pdbs:
        sys.exit(f"ERROR: no rank_001 PDB matching *{design}*rank_001*.pdb under {pred_dir}")
    pdb = pdbs[0]
    jsons = sorted(pred_dir.glob(f"*{design}*scores_rank_001*.json"))
    if not jsons:  # some ColabFold versions name it differently
        jsons = sorted(pred_dir.glob(f"*{design}*rank_001*.json"))
    if not jsons:
        sys.exit(f"ERROR: no rank_001 scores JSON for {design} under {pred_dir}")
    return pdb, jsons[0]


def _unwrap(d):
    return d[0] if isinstance(d, list) else d


def load_pae(json_path: Path) -> np.ndarray:
    d = _unwrap(json.loads(Path(json_path).read_text()))
    for key in ("pae", "predicted_aligned_error"):
        if key in d:
            return np.asarray(d[key], dtype=float)
    sys.exit(f"ERROR: no PAE key in {json_path} (keys: {list(d)[:8]})")


def load_plddt(json_path: Path):
    d = _unwrap(json.loads(Path(json_path).read_text()))
    return np.asarray(d["plddt"], dtype=float) if "plddt" in d else None


def ordered_residues(pdb_path: Path):
    """Flat residue list in PDB document order (== AF2 FASTA/PAE order)."""
    structure = PDBParser(QUIET=True).get_structure("m", str(pdb_path))
    model = next(structure.get_models())
    res = []
    for chain in model:
        for r in chain:
            if r.id[0] != " ":  # skip HETATM / water
                continue
            heavy = np.array([a.coord for a in r if a.element != "H"], dtype=float)
            try:
                aa1 = seq1(r.resname)
            except Exception:
                aa1 = "X"
            res.append(
                dict(chain=chain.id, resseq=r.id[1], aa1=aa1, heavy=heavy)
            )
    return res


def min_heavy_dist(a: np.ndarray, b: np.ndarray) -> float:
    if a.size == 0 or b.size == 0:
        return float("inf")
    diff = a[:, None, :] - b[None, :, :]
    return float(np.sqrt((diff * diff).sum(-1)).min())


def band(label, value, lo, hi):
    if value == float("inf"):
        return f"{label}: inf (NO contact <= 8 A)"
    inside = lo <= value <= hi
    return f"{label}: {value:.2f}  [{'in' if inside else 'OUT of'} positive band {lo}-{hi}]"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred-dir", default=DEFAULT_PRED_DIR)
    ap.add_argument("--design", default=DEFAULT_DESIGN)
    ap.add_argument("--binder-chain", default="D")
    ap.add_argument("--mhc-chains", default="A,B")
    args = ap.parse_args()

    pred_dir = Path(args.pred_dir)
    pdb, scores = find_files(pred_dir, args.design)
    print(f"PDB:  {pdb}")
    print(f"PAE:  {scores}\n")

    res = ordered_residues(pdb)
    pae = load_pae(scores)
    plddt = load_plddt(scores)

    if len(res) != pae.shape[0]:
        sys.exit(
            f"ERROR: residue count {len(res)} != PAE dim {pae.shape[0]}. "
            "Chain order / HETATM handling mismatch -- inspect before trusting output."
        )

    # --- chain identification (loud + verifiable, like the cross-check thread) ---
    chains = {}
    for i, r in enumerate(res):
        chains.setdefault(r["chain"], []).append(i)
    print("chains in model (id: len, seq-head):")
    pep_chain = None
    for cid, idxs in chains.items():
        seq = "".join(res[i]["aa1"] for i in idxs)
        tag = ""
        if seq == PEPTIDE_SEQ:
            pep_chain = cid
            tag = "   <- PEPTIDE (RMFPNAPYL)"
        elif cid == args.binder_chain:
            tag = "   <- BINDER"
        print(f"  {cid}: len={len(idxs):<4} {seq[:24]}{'...' if len(seq) > 24 else '':<3}{tag}")
    if pep_chain is None:
        sys.exit("ERROR: no chain matched RMFPNAPYL exactly. Inspect chain seqs above.")
    mhc_chains = [c.strip() for c in args.mhc_chains.split(",")]
    binder_chain = args.binder_chain
    print()

    binder_idx = chains[binder_chain]
    pep_idx = chains[pep_chain]
    mhc_idx = [i for c in mhc_chains for i in chains.get(c, [])]
    target_idx = pep_idx + mhc_idx

    # --- distance matrix: binder x target (heavy-atom min per residue pair) ---
    D = np.full((len(binder_idx), len(target_idx)), np.inf)
    for bi, b in enumerate(binder_idx):
        hb = res[b]["heavy"]
        for ti, t in enumerate(target_idx):
            D[bi, ti] = min_heavy_dist(hb, res[t]["heavy"])

    pep_cols = list(range(len(pep_idx)))
    mhc_cols = list(range(len(pep_idx), len(target_idx)))

    def iface(cols):
        scores_, contact_res = [], set()
        for bi, b in enumerate(binder_idx):
            for c in cols:
                if D[bi, c] <= CONTACT_A:
                    t = target_idx[c]
                    scores_.append(min(pae[b, t], pae[t, b]))
                    contact_res.add(b)
                    contact_res.add(t)
        val = float(np.mean(scores_)) if scores_ else float("inf")
        return val, contact_res

    iface_pep, _ = iface(pep_cols)
    iface_mhc, _ = iface(mhc_cols)
    iface_tot, contact_tot = iface(pep_cols + mhc_cols)

    # ipLDDT over interface-contact residues (tot)
    if plddt is not None and contact_tot:
        iplddt = float(np.mean([plddt[i] for i in contact_tot]))
    else:
        iplddt = float("nan")

    # --- per-peptide-position table ---
    print(" P aa role                  min_dist closest_binder_res PAE_to_binder  contact")
    print("-" * 82)
    exposed_contacted, anchor_contacted = [], []
    for p_local, b_pep in enumerate(pep_idx, start=1):
        dists = [(min_heavy_dist(res[b_pep]["heavy"], res[b]["heavy"]), res[b]["resseq"], b)
                 for b in binder_idx]
        dmin, closest_resseq, closest_flat = min(dists, key=lambda x: x[0])
        pae_pb = min(pae[b_pep, closest_flat], pae[closest_flat, b_pep])
        tier = "TIGHT" if dmin <= TIGHT_A else ("contact" if dmin <= CONTACT_A else "")
        aa = res[b_pep]["aa1"]
        role = PEP_ROLE.get(p_local, "")
        print(f"{p_local:2d}  {aa} {role:<22} {dmin:6.2f}  {closest_resseq:>16}  {pae_pb:11.2f}  {tier:>7}")
        if dmin <= CONTACT_A:
            if p_local in EXPOSED_POS:
                exposed_contacted.append(f"P{p_local}{aa}")
            if p_local in ANCHOR_POS:
                anchor_contacted.append(f"P{p_local}{aa}")

    binder_within5 = sorted({
        res[b]["resseq"]
        for bi, b in enumerate(binder_idx)
        for c in pep_cols
        if D[bi, c] <= TIGHT_A
    })

    # --- verdict ---
    print("\nVERDICT")
    print(f"  exposed/specificity residues contacted (<=8A): {exposed_contacted}")
    print(f"  anchor residues contacted (<=8A):              {anchor_contacted}")
    print(f"  binder residues within 5A of peptide:          {binder_within5}")
    print("\nINTERFACE-8A METRICS (comparable to controls + 3010)")
    print(f"  {band('iface_tot', iface_tot, *POS_TOT)}")
    print(f"  {band('iface_pep', iface_pep, *POS_PEP)}")
    print(f"  iface_mhc: {iface_mhc:.2f}" if iface_mhc != float('inf') else "  iface_mhc: inf")
    print(f"  ipLDDT (interface): {iplddt:.1f}")
    print("\n  reference  3010: tot 2.98 / pep 2.10 / mhc 3.35  (the peptide reader)")
    print("  reference  3054: tot 1.61 / pep inf              (framework champion)")
    pep_blind = iface_pep == float("inf") or iface_pep > POS_PEP[1] * 2
    print(
        "\n  -> 2079 reads the PEPTIDE"
        if not pep_blind
        else "\n  -> 2079 is PEPTIDE-BLIND (framework-only), like the c03 framework binders"
    )


if __name__ == "__main__":
    main()
