"""Unit tests for workflow.scripts.splice_binder."""

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


def _write_pmhc(path: Path) -> None:
    atoms = [
        ("A", 65, (0.0, 0.0, 0.0)),
        ("A", 66, (3.8, 0.0, 0.0)),
        ("A", 150, (0.0, 5.0, 0.0)),
        ("B", 1, (20.0, 0.0, 0.0)),
        ("B", 2, (20.0, 3.8, 0.0)),
        ("C", 1, (0.0, 0.0, 5.0)),
        ("C", 4, (3.0, 3.0, 5.0)),
        ("C", 9, (4.0, 9.0, 5.0)),
    ]
    _write_pdb(path, atoms)


def _write_stage1_binder(path: Path, resnums: list[int], chain: str = "A") -> None:
    atoms = [(chain, n, (10.0 + n * 0.1, 5.0, 0.0)) for n in resnums]
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
    _write_stage1_binder(s1, [10, 11, 12, 13, 14])
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
    _write_stage1_binder(s1, [42, 43, 44, 45])
    splice_binder_onto_pmhc(s1, pmhc, out)

    chains = _parse_chains(out)
    for cid, resnums in chains.items():
        assert resnums[0] == 1, f"chain {cid} should start at 1, got {resnums[0]}"
        assert resnums == list(range(1, len(resnums) + 1)), f"chain {cid} not contiguous: {resnums}"


def test_splice_rejects_input_with_chain_D(tmp_path: Path) -> None:
    pmhc = tmp_path / "pmhc.pdb"
    s1 = tmp_path / "design.pdb"
    out = tmp_path / "spliced.pdb"
    _write_pmhc(pmhc)
    _write_stage1_binder(s1, [1, 2, 3], chain="D")
    with pytest.raises(ValueError, match="expected binder on chain 'A'"):
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
    _write_stage1_binder(s1, [1, 2, 3])
    with pytest.raises(ValueError, match="chain 'D' already present"):
        splice_binder_onto_pmhc(s1, pmhc, out)


def test_splice_rejects_stage1_with_multiple_chains(tmp_path: Path) -> None:
    pmhc = tmp_path / "pmhc.pdb"
    s1 = tmp_path / "design.pdb"
    out = tmp_path / "spliced.pdb"
    _write_pmhc(pmhc)
    atoms = [
        ("A", 1, (10.0, 5.0, 0.0)),
        ("B", 1, (11.0, 5.0, 0.0)),
    ]
    _write_pdb(s1, atoms)
    with pytest.raises(ValueError, match="expected exactly 1 chain"):
        splice_binder_onto_pmhc(s1, pmhc, out)
