"""Deterministic fixtures for the cycle 03 BAKER-truncation mock tests.

Produces:
- tests/fixtures/targets/3hpj_baker_truncated_mock.pdb -- the BAKER-format
  truncated reference: chain B = HLA groove stub (180 CA, resnum 1-180),
  chain C = peptide (9 CA, resnum 1-9). Matches the real prep_baker_target.py
  output layout (no chain A, no beta2m) at full residue count so
  align_baker_scaffolds.py emits B=180 / C=9 in mock just as on the pod.
- tests/fixtures/baker_truncated_mock/scaf{0,1}.pdb -- mini BAKER scaffolds:
  chain A = binder (~25 CA), chain B = HLA(180) + peptide(9) fused as 189
  continuous residues (resnum 1-189). The HLA substring is the reference HLA
  rigidly translated by +100 A in x, so align_baker_scaffolds.py recovers a
  non-trivial transform (RMSD ~ 0).

Run from repo root:
    uv run python tests/fixtures/baker_truncated_mock/_make_fixtures.py
"""

from __future__ import annotations

from pathlib import Path

FIXTURE_DIR = Path(__file__).resolve().parent
TRUNCATED_TARGET = FIXTURE_DIR.parent / "targets" / "3hpj_baker_truncated_mock.pdb"

HLA_LEN = 180
PEPTIDE_LEN = 9

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


def _hla_coords(i: int) -> tuple[float, float, float]:
    """Deterministic spread-out HLA CA position for residue index i (0-based)."""
    return (i * 3.8, (i % 7) * 1.3, (i % 5) * 0.9)


def _peptide_coords(i: int) -> tuple[float, float, float]:
    return (i * 1.1, 5.0 + i * 0.4, 5.0)


def write_truncated_target() -> None:
    hla = [("B", i + 1, _hla_coords(i)) for i in range(HLA_LEN)]
    peptide = [("C", i + 1, _peptide_coords(i)) for i in range(PEPTIDE_LEN)]
    _write_pdb(TRUNCATED_TARGET, hla + peptide, "MOCK BAKER-TRUNCATED TARGET (B=HLA180, C=PEP9)")


def write_scaffolds() -> None:
    """Scaffolds: chain A binder + chain B fused HLA(+100 A in x) + peptide."""
    for idx in range(2):
        n_binder = 25 + idx
        binder = [("A", i + 1, (1.5, i * 3.9 + idx, 2.5)) for i in range(n_binder)]
        # Fused chain B: HLA residues 1-180 (reference + 100 A in x) then
        # peptide residues 181-189 at arbitrary nearby coords.
        fused: list[Atom] = []
        for i in range(HLA_LEN):
            x, y, z = _hla_coords(i)
            fused.append(("B", i + 1, (x + 100.0, y, z)))
        for i in range(PEPTIDE_LEN):
            fused.append(("B", HLA_LEN + i + 1, (100.0 + i, 6.0, 6.0)))
        _write_pdb(
            FIXTURE_DIR / f"scaf{idx}.pdb", binder + fused, f"BAKER MOCK FUSED SCAFFOLD {idx}"
        )


def main() -> None:
    write_truncated_target()
    write_scaffolds()
    print(f"wrote {TRUNCATED_TARGET} and scaffolds under {FIXTURE_DIR}")


if __name__ == "__main__":
    main()
