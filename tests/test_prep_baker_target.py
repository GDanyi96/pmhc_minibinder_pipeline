"""Unit tests for workflow.scripts.prep_baker_target (cycle 03 sub-run A).

The real input (data/targets/3hpj_clean.pdb) is produced pod-side, so these
tests synthesize a full canonical clean target (chain A = HC 275, B = beta2m
100, C = peptide 9) and assert the BAKER-format truncation produces exactly
chain B = HLA[1:180] and chain C = peptide[1:9].
"""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio.PDB import PDBParser

from workflow.scripts import prep_baker_target

REPO_ROOT = Path(__file__).resolve().parent.parent
MOCK_FIXTURE = REPO_ROOT / "tests/fixtures/targets/3hpj_baker_truncated_mock.pdb"


def _atom(serial: int, chain: str, resnum: int, xyz: tuple[float, float, float]) -> str:
    x, y, z = xyz
    return (
        f"ATOM  {serial:>5d}  CA  ALA {chain}{resnum:>4d}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C"
    )


def _write_full_clean(path: Path, hc: int = 275, b2m: int = 100, pep: int = 9) -> None:
    lines = ["HEADER    SYNTH FULL CLEAN 3HPJ"]
    serial = 1
    for i in range(hc):
        lines.append(_atom(serial, "A", i + 1, (i * 3.8, 0.0, 0.0)))
        serial += 1
    for i in range(b2m):
        lines.append(_atom(serial, "B", i + 1, (i * 3.8, 50.0, 0.0)))
        serial += 1
    for i in range(pep):
        lines.append(_atom(serial, "C", i + 1, (i * 1.1, 5.0, 5.0)))
        serial += 1
    # A HETATM water inside the kept chain-A range must be dropped.
    lines.append("HETATM 9999  O   HOH A  50      10.000  10.000  10.000  1.00  0.00           O")
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def _ca_counts(pdb: Path) -> dict[str, int]:
    model = PDBParser(QUIET=True).get_structure("x", str(pdb))[0]
    return {
        c.id: sum(1 for r in c.get_residues() if r.id[0].strip() == "" and "CA" in r)
        for c in model.get_chains()
    }


def test_truncation_yields_180_b_and_9_c(tmp_path: Path) -> None:
    src = tmp_path / "3hpj_clean.pdb"
    _write_full_clean(src)
    out = tmp_path / "3hpj_baker_truncated.pdb"
    counts = prep_baker_target.truncate_target(src, out)
    assert counts == {"B": 180, "C": 9}
    assert _ca_counts(out) == {"B": 180, "C": 9}


def test_no_chain_a_or_beta2m_and_numbering_preserved(tmp_path: Path) -> None:
    src = tmp_path / "3hpj_clean.pdb"
    _write_full_clean(src)
    out = tmp_path / "trunc.pdb"
    prep_baker_target.truncate_target(src, out)
    model = PDBParser(QUIET=True).get_structure("x", str(out))[0]
    assert {c.id for c in model.get_chains()} == {"B", "C"}
    b_resnums = [r.id[1] for r in model["B"].get_residues() if r.id[0].strip() == ""]
    assert b_resnums == list(range(1, 181))  # original 1..180 preserved, alpha3 dropped
    # No HETATM survives.
    assert all(r.id[0].strip() == "" for c in model.get_chains() for r in c.get_residues())


def test_rejects_missing_chains(tmp_path: Path) -> None:
    src = tmp_path / "bad.pdb"
    src.write_text("HEADER  X\n" + _atom(1, "A", 1, (0, 0, 0)) + "\nEND\n")
    with pytest.raises(ValueError, match="absent"):
        prep_baker_target.truncate_target(src, tmp_path / "out.pdb")


def test_rejects_short_groove(tmp_path: Path) -> None:
    src = tmp_path / "short.pdb"
    _write_full_clean(src, hc=150)  # alpha1+alpha2 incomplete
    with pytest.raises(ValueError, match="HLA groove"):
        prep_baker_target.truncate_target(src, tmp_path / "out.pdb")


def test_mock_copies_fixture(tmp_path: Path) -> None:
    out = tmp_path / "trunc.pdb"
    rc = prep_baker_target.main(["--mock", "--out", str(out)])
    assert rc == 0
    assert _ca_counts(out) == {"B": 180, "C": 9}


def test_mock_fixture_has_expected_layout() -> None:
    assert _ca_counts(MOCK_FIXTURE) == {"B": 180, "C": 9}
