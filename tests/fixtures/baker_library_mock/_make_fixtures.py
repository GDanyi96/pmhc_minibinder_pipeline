"""Deterministic fixture generator for the Cycle 03 prep mock tests.

Produces:
- tests/fixtures/baker_library_mock/scaf{0,1,2}.pdb -- 3 synthetic mini
  BAKER-style scaffolds. Each has chain A = binder (~20 CA), chain B =
  truncated MHC target stub (4 CA matching the reference MHC chain A but
  rigidly translated by +100 A in x, so align_scaffolds.py recovers a
  non-trivial transform). Mirrors BAKER's scaf*.pdb layout (binder + target
  context) at toy scale.
- tests/fixtures/design_2079_mock_seed.pdb -- the sub-run B seed: a full
  4-chain complex (A/B/C copied from the reference mock_clean.pdb + chain D
  = a ~20-residue 4-helix-bundle stand-in for the cycle-02 hero), stitched
  exactly as workflow/scripts/setup_cycle03_inputs.py builds the real seed.

Reference frame is tests/fixtures/stage1/mock_clean.pdb (chains A/B/C,
hotspot cluster). Binder backbones are placed through the hotspot-cluster
centroid so the merged Stage 1 designs pass the geometry gate.

Run from repo root:
    uv run python tests/fixtures/baker_library_mock/_make_fixtures.py
"""

from __future__ import annotations

from pathlib import Path

FIXTURE_DIR = Path(__file__).resolve().parent
SEED_OUT = FIXTURE_DIR.parent / "design_2079_mock_seed.pdb"

# A/B/C motif atoms -- identical to tests/fixtures/stage1/mock_clean.pdb so
# the spliced/merged complexes share the reference target frame.
MHC_A = [
    ("A", 65, (0.0, 0.0, 0.0)),
    ("A", 66, (3.8, 0.0, 0.0)),
    ("A", 150, (0.0, 5.0, 0.0)),
    ("A", 155, (0.0, 8.8, 0.0)),
]
B2M_B = [
    ("B", 1, (20.0, 0.0, 0.0)),
    ("B", 2, (20.0, 3.8, 0.0)),
]
PEPTIDE_C = [
    ("C", 1, (0.0, 0.0, 5.0)),
    ("C", 2, (1.0, 1.0, 5.0)),
    ("C", 3, (2.0, 2.0, 5.0)),
    ("C", 4, (3.0, 3.0, 5.0)),
    ("C", 5, (1.5, 4.0, 5.0)),
    ("C", 6, (0.0, 5.0, 5.0)),
    ("C", 7, (1.5, 6.5, 5.0)),
    ("C", 8, (3.0, 8.0, 5.0)),
    ("C", 9, (4.0, 9.0, 5.0)),
]

Atom = tuple[str, int, tuple[float, float, float]]


def _atom_line(serial: int, chain: str, resnum: int, xyz: tuple[float, float, float]) -> str:
    x, y, z = xyz
    return (
        f"ATOM  {serial:>5d}  CA  ALA {chain}{resnum:>4d}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C"
    )


def _write_pdb(path: Path, atoms: list[Atom], header: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"HEADER    {header}"]
    for serial, (chain, resnum, xyz) in enumerate(atoms, start=1):
        lines.append(_atom_line(serial, chain, resnum, xyz))
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def _binder_backbone(chain: str, n: int) -> list[Atom]:
    """Straight backbone through the hotspot-cluster centroid (passes the
    Stage 1 geometry gate: >=3 CA within 10 A of a hotspot)."""
    atoms: list[Atom] = []
    for i in range(n):
        x = 1.5
        y = i * 3.9 - (n - 1) * 3.9 / 2 + 3.5
        z = 2.5
        atoms.append((chain, i + 1, (x, y, z)))
    return atoms


def write_scaffolds() -> None:
    """3 scaffolds: chain A binder + chain B truncated target (+100 A in x)."""
    for idx in range(3):
        binder = _binder_backbone("A", 18 + idx)  # arbitrary scaffold binder
        target_stub: list[Atom] = [
            ("B", resnum, (x + 100.0, y, z)) for (_c, resnum, (x, y, z)) in MHC_A
        ]
        _write_pdb(
            FIXTURE_DIR / f"scaf{idx}.pdb",
            binder + target_stub,
            f"BAKER MOCK SCAFFOLD {idx}",
        )


def write_seed_complex() -> None:
    """design_2079 mock seed: A/B/C reference motif + chain D 4HB binder."""
    binder_d = _binder_backbone("D", 20)
    _write_pdb(
        SEED_OUT,
        MHC_A + B2M_B + PEPTIDE_C + binder_d,
        "DESIGN_2079 MOCK SEED (STITCHED COMPLEX)",
    )


def main() -> None:
    write_scaffolds()
    write_seed_complex()
    print(f"wrote scaffolds under {FIXTURE_DIR} and seed at {SEED_OUT}")


if __name__ == "__main__":
    main()
