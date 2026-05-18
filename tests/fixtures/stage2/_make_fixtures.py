"""Generate deterministic ColabFold mock fixtures for the 5 Stage 2 controls.

Run from the repo root:
    python tests/fixtures/stage2/_make_fixtures.py

For each control the JSON encodes a 2-D PAE matrix sized to its chain
concatenation length (A=HC, B=beta2m, C=peptide, D=binder). Within a chain
the PAE is 0.5; the binder<->pMHC submatrix carries a tuned value so iPAE
comes out close to the target.

Targets match cycle-1 pod observations (avoid mock-vs-real divergence):
P1 iPAE 5.0 / ipLDDT 95.0  (real: 4.88 / 94.5)
P2 iPAE 6.5 / ipLDDT 95.0  (real: 4.56 / 95.9; Jenkins published ~6.5)
P3 iPAE 4.5 / ipLDDT 95.0  (real: 4.54 / 94.2; ~ P1 -- AF2 blind spot)
N1 iPAE 25.0 / ipLDDT 35.0 (real: 24.72 / 34.5)
N2 iPAE 27.0 / ipLDDT 35.0 (real: 25.70 / 31.4)

These values exercise the commit-5 halt rule (positives <-> negatives
iPAE gap > 10 Ang and ipLDDT gap > 30) at realistic separation, so the
mock catches regressions that would matter on real data.

JSON key is "pae" (ColabFold >=1.6.0). compute_metrics.load_pae also
accepts the legacy "predicted_aligned_error" key. ptm/iptm are written
but unused by our code (pitfall #2).

Reproducible: no randomness, idempotent on rerun.
"""

from __future__ import annotations

import json
from pathlib import Path

FIXTURE_DIR = Path(__file__).resolve().parent / "controls_colabfold"

CONTROLS = {
    "P1": {"hc": 275, "b2m": 100, "pep": 9, "binder": 85, "ipae_target": 5.0, "plddt": 95.0},
    "P2": {"hc": 276, "b2m": 100, "pep": 9, "binder": 145, "ipae_target": 6.5, "plddt": 95.0},
    "P3": {"hc": 275, "b2m": 100, "pep": 10, "binder": 85, "ipae_target": 4.5, "plddt": 95.0},
    "N1": {"hc": 275, "b2m": 100, "pep": 9, "binder": 85, "ipae_target": 25.0, "plddt": 35.0},
    "N2": {"hc": 275, "b2m": 100, "pep": 9, "binder": 85, "ipae_target": 27.0, "plddt": 35.0},
}

SCORES_TEMPLATE = "{cid}_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json"
PDB_TEMPLATE = "{cid}_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb"


def build_pae(
    n: int, binder_slice: slice, mhc_slice: slice, cross_value: float
) -> list[list[float]]:
    matrix = [[0.5] * n for _ in range(n)]
    for i in range(binder_slice.start, binder_slice.stop):
        for j in range(mhc_slice.start, mhc_slice.stop):
            matrix[i][j] = cross_value
            matrix[j][i] = cross_value
    return matrix


def minimal_pdb(chain_lengths: dict[str, int]) -> str:
    lines: list[str] = []
    atom_serial = 1
    res_serial = 1
    # One CA per residue, evenly spaced along x within a chain, chains offset on y.
    for ci, (chain_id, length) in enumerate(chain_lengths.items()):
        for r in range(length):
            x = float(r)
            y = float(ci * 10)
            z = 0.0
            lines.append(
                f"ATOM  {atom_serial:>5d}  CA  GLY {chain_id}{res_serial:>4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 95.00           C"
            )
            atom_serial += 1
            res_serial += 1
        lines.append(f"TER   {atom_serial:>5d}      GLY {chain_id}{res_serial - 1:>4d}")
        atom_serial += 1
        res_serial = 1
    lines.append("END")
    return "\n".join(lines) + "\n"


def write_fixtures() -> None:
    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    for cid, lens in CONTROLS.items():
        n = lens["hc"] + lens["b2m"] + lens["pep"] + lens["binder"]
        a = slice(0, lens["hc"])
        b = slice(a.stop, a.stop + lens["b2m"])
        c = slice(b.stop, b.stop + lens["pep"])
        d = slice(c.stop, c.stop + lens["binder"])
        mhc_total = (a.stop - a.start) + (b.stop - b.start) + (c.stop - c.start)
        binder_total = d.stop - d.start
        # Cross-chain block (binder<->pMHC, both directions) dominates the iPAE
        # mean; within-chain entries are 0.5 and not summed by compute_ipae.
        # Therefore iPAE == cross_value exactly when the block is uniform.
        cross = float(lens["ipae_target"])
        pae = [[0.5] * n for _ in range(n)]
        for i in range(d.start, d.stop):
            for j in range(a.start, c.stop):
                pae[i][j] = cross
                pae[j][i] = cross

        scores = {
            # ColabFold >=1.6.0 emits the matrix under "pae".
            "pae": pae,
            "plddt": [float(lens["plddt"])] * n,
            "ptm": 0.85,
            "iptm": 0.75,
            "max_pae": cross,
        }
        json_path = FIXTURE_DIR / SCORES_TEMPLATE.format(cid=cid)
        json_path.write_text(json.dumps(scores))

        chain_lengths = {
            "A": lens["hc"],
            "B": lens["b2m"],
            "C": lens["pep"],
            "D": lens["binder"],
        }
        pdb_path = FIXTURE_DIR / PDB_TEMPLATE.format(cid=cid)
        pdb_path.write_text(minimal_pdb(chain_lengths))
        # Echo so reruns are auditable.
        print(
            f"{cid}: N={n} cross={cross} binder_atoms={binder_total} mhc_atoms={mhc_total} "
            f"-> {json_path.name}, {pdb_path.name}"
        )


if __name__ == "__main__":
    write_fixtures()
