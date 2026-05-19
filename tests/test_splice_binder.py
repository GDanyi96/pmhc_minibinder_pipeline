"""Unit tests for workflow.scripts.splice_binder.

Stage 1 (RFdiffusion) emits a 4-chain complex (A=HC, B=beta2m, C=peptide,
D=binder); the splicer extracts chain D and combines it with A/B/C from
the canonical cleaned pMHC. See docs/known_traps.md trap #29.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio.PDB import PDBParser

from workflow.scripts.splice_binder import splice_binder_onto_pmhc


def _atom_line(serial: int, chain: str, resnum: int, xyz: tuple[float, float, float]) -> str:
    x, y, z = xyz
    return (
        f"ATOM  {serial:>5d}  CA  ALA {chain}{resnum:>4d}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C"
    )


def _write_pdb(path: Path, atoms: list[tuple[str, int, tuple[float, float, float]]]) -> None:
    lines = ["HEADER    SPLICE TEST FIXTURE"]
    for i, (chain, resnum, xyz) in enumerate(atoms, start=1):
        lines.append(_atom_line(i, chain, resnum, xyz))
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def _pmhc_atoms() -> list[tuple[str, int, tuple[float, float, float]]]:
    return [
        ("A", 65, (0.0, 0.0, 0.0)),
        ("A", 66, (3.8, 0.0, 0.0)),
        ("A", 150, (0.0, 5.0, 0.0)),
        ("B", 1, (20.0, 0.0, 0.0)),
        ("B", 2, (20.0, 3.8, 0.0)),
        ("C", 1, (0.0, 0.0, 5.0)),
        ("C", 4, (3.0, 3.0, 5.0)),
        ("C", 9, (4.0, 9.0, 5.0)),
    ]


def _write_pmhc(path: Path) -> None:
    _write_pdb(path, _pmhc_atoms())


def _write_stage1_4chain(
    path: Path,
    binder_resnums: list[int],
    binder_chain: str = "D",
) -> None:
    """Write a 4-chain Stage 1 PDB. A/B/C mirror the cleaned pMHC; the
    binder is written to ``binder_chain`` (default D)."""
    atoms = list(_pmhc_atoms())
    for n in binder_resnums:
        atoms.append((binder_chain, n, (10.0 + n * 0.1, 5.0, 0.0)))
    _write_pdb(path, atoms)


def _parse_chains(pdb: Path) -> dict[str, list[int]]:
    parser = PDBParser(QUIET=True)
    model = parser.get_structure(pdb.stem, str(pdb))[0]
    return {c.id: [r.id[1] for r in c.get_residues()] for c in model.get_chains()}


def test_splice_chain_renaming(tmp_path: Path) -> None:
    pmhc = tmp_path / "pmhc.pdb"
    s1 = tmp_path / "design.pdb"
    out = tmp_path / "spliced.pdb"
    _write_pmhc(pmhc)
    _write_stage1_4chain(s1, [10, 11, 12, 13, 14])
    splice_binder_onto_pmhc(s1, pmhc, out)

    chains = _parse_chains(out)
    assert set(chains.keys()) == {"A", "B", "C", "D"}, chains
    assert len(chains["A"]) == 3
    assert len(chains["B"]) == 2
    assert len(chains["C"]) == 3
    assert len(chains["D"]) == 5


def test_splice_residue_renumbering(tmp_path: Path) -> None:
    pmhc = tmp_path / "pmhc.pdb"
    s1 = tmp_path / "design.pdb"
    out = tmp_path / "spliced.pdb"
    _write_pmhc(pmhc)
    _write_stage1_4chain(s1, [42, 43, 44, 45])
    splice_binder_onto_pmhc(s1, pmhc, out)

    chains = _parse_chains(out)
    for cid, resnums in chains.items():
        assert resnums[0] == 1, f"chain {cid} should start at 1, got {resnums[0]}"
        assert resnums == list(range(1, len(resnums) + 1)), f"chain {cid} not contiguous: {resnums}"


def test_splice_4chain_stage1_input_uses_canonical_abc(tmp_path: Path) -> None:
    """The output A/B/C must come from the cleaned pMHC, not the Stage 1
    PDB. We verify by writing a Stage 1 PDB whose A/B/C residue counts
    differ from the cleaned pMHC's, and confirming the output A/B/C
    matches the cleaned pMHC counts.
    """
    pmhc = tmp_path / "pmhc.pdb"
    s1 = tmp_path / "design.pdb"
    out = tmp_path / "spliced.pdb"
    _write_pmhc(pmhc)  # A=3, B=2, C=3 residues

    # 4-chain Stage 1 PDB with deliberately-different A/B/C counts (an
    # abbreviated stand-in for the cycle 02 contig A=275, B=100, C=9,
    # D=85). All four chains are present so the unexpected-chain guard
    # is not tripped.
    s1_atoms: list[tuple[str, int, tuple[float, float, float]]] = []
    for resnum in range(1, 11):  # A=10
        s1_atoms.append(("A", resnum, (resnum * 0.5, 0.0, 0.0)))
    for resnum in range(1, 6):  # B=5
        s1_atoms.append(("B", resnum, (20.0, resnum * 0.5, 0.0)))
    for resnum in range(1, 5):  # C=4
        s1_atoms.append(("C", resnum, (0.0, 0.0, resnum * 0.5)))
    binder_resnums = list(range(1, 8))  # D=7
    for resnum in binder_resnums:
        s1_atoms.append(("D", resnum, (10.0 + resnum * 0.1, 5.0, 0.0)))
    _write_pdb(s1, s1_atoms)

    splice_binder_onto_pmhc(s1, pmhc, out)
    chains = _parse_chains(out)
    assert set(chains.keys()) == {"A", "B", "C", "D"}
    # A/B/C from cleaned pMHC (3/2/3), NOT from Stage 1 (10/5/4)
    assert len(chains["A"]) == 3
    assert len(chains["B"]) == 2
    assert len(chains["C"]) == 3
    # D from Stage 1
    assert len(chains["D"]) == len(binder_resnums)


def test_splice_rejects_input_without_chain_D(tmp_path: Path) -> None:
    """Stage 1 PDB without chain D (e.g. the cycle 01 single-chain
    layout) must fail loudly rather than silently re-mapping a
    different chain."""
    pmhc = tmp_path / "pmhc.pdb"
    s1 = tmp_path / "design.pdb"
    out = tmp_path / "spliced.pdb"
    _write_pmhc(pmhc)
    # Single chain "A", no D.
    _write_pdb(s1, [("A", n, (10.0 + n * 0.1, 5.0, 0.0)) for n in (1, 2, 3)])
    with pytest.raises(ValueError, match="expected binder on chain 'D'"):
        splice_binder_onto_pmhc(s1, pmhc, out)


def test_splice_rejects_unexpected_chain(tmp_path: Path) -> None:
    """Stage 1 PDB with a chain outside {A,B,C,D} (e.g. an extra "E"
    from a future contig change) must fail loudly. Same defensive
    pattern as ``_binder_ca_coords`` (Trap #28)."""
    pmhc = tmp_path / "pmhc.pdb"
    s1 = tmp_path / "design.pdb"
    out = tmp_path / "spliced.pdb"
    _write_pmhc(pmhc)
    atoms = list(_pmhc_atoms())
    atoms.append(("D", 1, (10.0, 5.0, 0.0)))
    atoms.append(("E", 1, (15.0, 5.0, 0.0)))
    _write_pdb(s1, atoms)
    with pytest.raises(ValueError, match="unexpected chain"):
        splice_binder_onto_pmhc(s1, pmhc, out)


def test_splice_rejects_pmhc_with_chain_D(tmp_path: Path) -> None:
    pmhc = tmp_path / "pmhc.pdb"
    s1 = tmp_path / "design.pdb"
    out = tmp_path / "spliced.pdb"
    atoms = [
        ("A", 65, (0.0, 0.0, 0.0)),
        ("B", 1, (20.0, 0.0, 0.0)),
        ("C", 1, (0.0, 0.0, 5.0)),
        ("D", 1, (10.0, 5.0, 0.0)),
    ]
    _write_pdb(pmhc, atoms)
    _write_stage1_4chain(s1, [1, 2, 3])
    with pytest.raises(ValueError, match="chain 'D' already present"):
        splice_binder_onto_pmhc(s1, pmhc, out)
