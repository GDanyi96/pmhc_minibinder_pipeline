"""Generate deterministic ColabFold mock fixtures for the 5 Stage 2 controls.

Run from the repo root:
    python tests/fixtures/stage2/_make_fixtures.py

For each control the JSON encodes a 2-D PAE matrix sized to its chain
concatenation length (A=HC, B=beta2m, C=peptide, D=binder). Within a chain
the PAE is 0.5; the binder<->pMHC submatrix carries a tuned value so that
iPAE comes out close to the target (P1=5.0, P2=6.5, P3=9.0, N1=17.0,
N2=22.0). pLDDT is uniformly 95.0. ptm/iptm are written but unused by our
code (pitfall #2).

Reproducible: no randomness, idempotent on rerun.
"""

from __future__ import annotations

import json
from pathlib import Path

FIXTURE_DIR = Path(__file__).resolve().parent / "controls_colabfold"

CONTROLS = {
    "P1": {"hc": 275, "b2m": 100, "pep": 9, "binder": 85, "ipae_target": 5.0},
    "P2": {"hc": 276, "b2m": 100, "pep": 9, "binder": 145, "ipae_target": 6.5},
    "P3": {"hc": 275, "b2m": 100, "pep": 10, "binder": 85, "ipae_target": 9.0},
    "N1": {"hc": 275, "b2m": 100, "pep": 9, "binder": 85, "ipae_target": 17.0},
    "N2": {"hc": 275, "b2m": 100, "pep": 9, "binder": 85, "ipae_target": 22.0},
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
            "predicted_aligned_error": pae,
            "plddt": [95.0] * n,
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
