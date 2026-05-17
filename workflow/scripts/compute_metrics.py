# mypy: disable-error-code="no-untyped-call,attr-defined,no-any-return,misc"
"""Stage 2c — iPAE, ipLDDT, and BSA computation from ColabFold outputs.

iPAE definition (HADRUP_JENKINS_2025 Fig 1B): symmetric mean PAE between
chain D (binder) and chains A+B+C (HC + beta2m + peptide).

ipLDDT definition: mean pLDDT over chain D residues.

BSA: SASA(pMHC) + SASA(binder) - SASA(complex) via freesasa default
Lee-Richards. Splits the complex PDB into A+B+C and D-only temporary PDBs
written via Biopython PDBIO so freesasa.Structure.from_file consumes only
ATOM records (pitfall #8).

Chain order A=HC, B=beta2m, C=peptide, D=binder. Boundaries are computed
from the actual PDB at runtime, never hard-coded (pitfall #4). The JSON
key `iptm` is never used as an iPAE proxy (pitfall #2).
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from glob import glob
from pathlib import Path
from typing import Any

import freesasa
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure

CHAIN_ORDER = ("A", "B", "C", "D")
DEFAULT_SCORES_GLOB = "*_scores_rank_001*.json"
DEFAULT_PDB_GLOB = "*_unrelaxed_rank_001*.pdb"


def load_pae(scores_json: Path) -> list[list[float]]:
    payload = json.loads(scores_json.read_text())
    pae = payload.get("predicted_aligned_error")
    if pae is None:
        raise KeyError(f"{scores_json}: missing 'predicted_aligned_error'")
    return pae


def load_plddt(scores_json: Path) -> list[float]:
    payload = json.loads(scores_json.read_text())
    plddt = payload.get("plddt")
    if plddt is None:
        raise KeyError(f"{scores_json}: missing 'plddt'")
    return plddt


def chain_lengths(pdb_path: Path) -> dict[str, int]:
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = structure[0]
    lengths: dict[str, int] = {}
    for chain in model.get_chains():
        if chain.id in CHAIN_ORDER:
            lengths[chain.id] = sum(1 for r in chain.get_residues() if r.id[0].strip() == "")
    return lengths


def chain_slices(lengths: dict[str, int]) -> dict[str, slice]:
    slices: dict[str, slice] = {}
    cursor = 0
    for chain_id in CHAIN_ORDER:
        size = lengths.get(chain_id, 0)
        slices[chain_id] = slice(cursor, cursor + size)
        cursor += size
    return slices


def compute_ipae(pae: list[list[float]], slices: dict[str, slice]) -> float:
    d_range = range(slices["D"].start, slices["D"].stop)
    mhc_ranges: list[range] = [range(slices[c].start, slices[c].stop) for c in ("A", "B", "C")]
    total = 0.0
    count = 0
    for d in d_range:
        for r in mhc_ranges:
            for j in r:
                total += pae[d][j] + pae[j][d]
                count += 2
    if count == 0:
        raise ValueError("compute_ipae: empty cross-chain submatrix (check chain ordering)")
    return total / count


def compute_iplddt(plddt: list[float], slices: dict[str, slice]) -> float:
    binder = plddt[slices["D"]]
    if not binder:
        raise ValueError("compute_iplddt: chain D is empty (check chain ordering)")
    return sum(binder) / len(binder)


class _DropChain(Select):
    def __init__(self, drop: str) -> None:
        self._drop = drop

    def accept_chain(self, chain: Chain) -> bool:
        return chain.id != self._drop


class _KeepChain(Select):
    def __init__(self, keep: str) -> None:
        self._keep = keep

    def accept_chain(self, chain: Chain) -> bool:
        return chain.id == self._keep


def _write_subset(structure: Structure, selector: Select, path: Path) -> None:
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(path), select=selector)


def _sasa(pdb_path: Path) -> float:
    structure = freesasa.Structure(str(pdb_path))
    result = freesasa.calc(structure)
    return float(result.totalArea())


def compute_bsa(complex_pdb: Path) -> float:
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(complex_pdb.stem, str(complex_pdb))
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        pmhc_pdb = tmp_dir / "pmhc.pdb"
        binder_pdb = tmp_dir / "binder.pdb"
        _write_subset(structure, _DropChain("D"), pmhc_pdb)
        _write_subset(structure, _KeepChain("D"), binder_pdb)
        sasa_complex = _sasa(complex_pdb)
        sasa_pmhc = _sasa(pmhc_pdb)
        sasa_binder = _sasa(binder_pdb)
    return sasa_pmhc + sasa_binder - sasa_complex


def find_artifact(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if not matches:
        matches = sorted(Path(p) for p in glob(str(directory / pattern)))
    if not matches:
        raise FileNotFoundError(f"no file matching {pattern!r} under {directory}")
    return matches[0]


def metrics_for_control(colabfold_dir: Path, backbone_pdb: Path | None) -> dict[str, Any]:
    scores_json = find_artifact(colabfold_dir, DEFAULT_SCORES_GLOB)
    rank_pdb = find_artifact(colabfold_dir, DEFAULT_PDB_GLOB)
    # Chain lengths come from the predicted complex PDB itself — the cleaned
    # crystal structure (when provided) lacks chain D, so we always slice from
    # the ColabFold rank_001 PDB, which has the full A+B+C+D layout.
    lengths = chain_lengths(rank_pdb)
    slices = chain_slices(lengths)
    pae = load_pae(scores_json)
    plddt = load_plddt(scores_json)
    ipae = compute_ipae(pae, slices)
    iplddt = compute_iplddt(plddt, slices)
    bsa = compute_bsa(rank_pdb)
    return {
        "ipae": ipae,
        "iplddt": iplddt,
        "bsa": bsa,
        "rank_001_pdb": str(rank_pdb),
        "scores_json": str(scores_json),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--colabfold-dir", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--mock", action="store_true")
    # Back-compat with original stub signature.
    parser.add_argument("--config")
    args = parser.parse_args(argv)

    colabfold_dir = args.colabfold_dir or (Path(args.config).parent if args.config else None)
    out_path = args.out
    if colabfold_dir is None or out_path is None:
        parser.error("must provide --colabfold-dir and --out (or legacy --config/--out)")

    if args.mock:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(
            json.dumps({"mock": True, "controls": {}, "designs": []}, indent=2) + "\n"
        )
        print(f"[mock] compute_metrics: wrote stub {out_path}")
        return 0

    result = metrics_for_control(colabfold_dir, None)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
