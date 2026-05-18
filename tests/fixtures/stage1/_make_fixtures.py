"""Deterministic fixture generator for Stage 1 mock tests.

Produces:
- mock_clean.pdb        synthetic pMHC cleaned PDB with chains A/B/C
                        containing all eight hotspot residues
                        (A65,A66,A150,A155 + C1,C4,C6,C8).
- mock_target.yaml      TargetManifest pointing at mock_clean.pdb.
- mock_design_NNNNN.pdb 10 binder backbones (8 geometry-pass / 2 fail).

Run from repo root:
    uv run python tests/fixtures/stage1/_make_fixtures.py
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

FIXTURE_DIR = Path(__file__).resolve().parent

HOTSPOTS_PEPTIDE = [
    ("C", 1, (0.0, 0.0, 5.0)),
    ("C", 4, (3.0, 3.0, 5.0)),
    ("C", 6, (0.0, 5.0, 5.0)),
    ("C", 8, (3.0, 8.0, 5.0)),
]
HOTSPOTS_MHC = [
    ("A", 65, (0.0, 0.0, 0.0)),
    ("A", 66, (3.8, 0.0, 0.0)),
    ("A", 150, (0.0, 5.0, 0.0)),
    ("A", 155, (0.0, 8.8, 0.0)),
]
PEPTIDE_FILLERS = [
    ("C", 2, (1.0, 1.0, 5.0)),
    ("C", 3, (2.0, 2.0, 5.0)),
    ("C", 5, (1.5, 4.0, 5.0)),
    ("C", 7, (1.5, 6.5, 5.0)),
    ("C", 9, (4.0, 9.0, 5.0)),
]
B2M_FILLERS = [
    ("B", 1, (20.0, 0.0, 0.0)),
    ("B", 2, (20.0, 3.8, 0.0)),
]

PASS_LENGTHS = [78, 82, 85, 87, 89, 92, 98, 104]
FAIL_LENGTHS = [50, 80]  # 50 = length oob; 80 = placed far from hotspots


def _atom_line(serial: int, chain: str, resnum: int, xyz: tuple[float, float, float]) -> str:
    x, y, z = xyz
    return (
        f"ATOM  {serial:>5d}  CA  ALA {chain}{resnum:>4d}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C"
    )


def _write_pdb(path: Path, atoms: list[tuple[str, int, tuple[float, float, float]]]) -> None:
    lines = ["HEADER    MOCK FIXTURE FOR STAGE 1"]
    for i, (chain, resnum, xyz) in enumerate(atoms, start=1):
        lines.append(_atom_line(i, chain, resnum, xyz))
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def write_mock_clean_pdb(out: Path) -> None:
    """Synthetic cleaned PDB: chain A (4 hotspot residues), chain B (2),
    chain C (9, including 4 hotspots + 5 fillers)."""
    atoms: list[tuple[str, int, tuple[float, float, float]]] = []
    atoms.extend(HOTSPOTS_MHC)
    atoms.extend(B2M_FILLERS)
    pep = sorted(HOTSPOTS_PEPTIDE + PEPTIDE_FILLERS, key=lambda a: a[1])
    atoms.extend(pep)
    _write_pdb(out, atoms)


def write_mock_target_yaml(out: Path, cleaned_pdb: Path) -> None:
    """Pre-built TargetManifest pointing at cleaned_pdb.

    cleaned_pdb is written as a bare filename relative to the manifest's
    own directory. Absolute paths in committed fixtures are CI-fragile
    (they bake in the generator's sandbox path and silently break when
    the repo is cloned elsewhere) -- see docs/known_traps.md trap #16.
    """
    primary = {
        "target_id": "wt1_a0201",
        "pdb_id": "3HPJ",
        "cleaned_pdb": cleaned_pdb.name,
        "peptide_sequence": "RMFPNAPYL",
        "chains": {
            "A": {"role": "heavy_chain", "length": 4},
            "B": {"role": "beta2m", "length": 2},
            "C": {"role": "peptide", "length": 9},
            "D": {"role": "binder", "length": None},
        },
        "hotspots": {
            "peptide": [
                {"token": "C1", "chain": "C", "resnum_clean": 1, "resnum_rcsb": 1},
                {"token": "C4", "chain": "C", "resnum_clean": 4, "resnum_rcsb": 4},
                {"token": "C6", "chain": "C", "resnum_clean": 6, "resnum_rcsb": 6},
                {"token": "C8", "chain": "C", "resnum_clean": 8, "resnum_rcsb": 8},
            ],
            "mhc": [
                {"token": "A65", "chain": "A", "resnum_clean": 65, "resnum_rcsb": 65},
                {"token": "A66", "chain": "A", "resnum_clean": 66, "resnum_rcsb": 66},
                {"token": "A150", "chain": "A", "resnum_clean": 150, "resnum_rcsb": 150},
                {"token": "A155", "chain": "A", "resnum_clean": 155, "resnum_rcsb": 155},
            ],
        },
        "binder": {"length_min": 70, "length_max": 110},
    }
    manifest: dict[str, Any] = {
        "cycle": "99",
        "mock": True,
        "primary": primary,
        # positive_control mirrors primary; Stage 1 only consumes .primary.
        "positive_control": primary,
    }
    out.write_text(yaml.safe_dump(manifest, sort_keys=False))


def _passing_backbone(n: int) -> list[tuple[str, int, tuple[float, float, float]]]:
    """Straight backbone passing through the hotspot cluster centroid."""
    atoms: list[tuple[str, int, tuple[float, float, float]]] = []
    for i in range(n):
        x = 1.5
        y = i * 3.9 - (n - 1) * 3.9 / 2 + 3.5
        z = 2.5
        atoms.append(("A", i + 1, (x, y, z)))
    return atoms


def _distant_backbone(n: int) -> list[tuple[str, int, tuple[float, float, float]]]:
    """Backbone placed 100 A away from the hotspot cluster — no contacts."""
    atoms: list[tuple[str, int, tuple[float, float, float]]] = []
    for i in range(n):
        atoms.append(("A", i + 1, (100.0, float(i) * 3.9, 100.0)))
    return atoms


def write_designs(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    idx = 0
    for n in PASS_LENGTHS:
        _write_pdb(out_dir / f"mock_design_{idx:05d}.pdb", _passing_backbone(n))
        idx += 1
    # 50-aa passes the contact test but fails length_out_of_range.
    _write_pdb(out_dir / f"mock_design_{idx:05d}.pdb", _passing_backbone(FAIL_LENGTHS[0]))
    idx += 1
    # 80-aa placed far from hotspots fails insufficient_hotspot_contacts.
    _write_pdb(out_dir / f"mock_design_{idx:05d}.pdb", _distant_backbone(FAIL_LENGTHS[1]))


def main() -> None:
    cleaned_pdb = FIXTURE_DIR / "mock_clean.pdb"
    write_mock_clean_pdb(cleaned_pdb)
    write_mock_target_yaml(FIXTURE_DIR / "mock_target.yaml", cleaned_pdb)
    write_designs(FIXTURE_DIR)
    print(f"wrote fixtures under {FIXTURE_DIR}")


if __name__ == "__main__":
    main()
