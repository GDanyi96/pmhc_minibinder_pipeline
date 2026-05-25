"""Unit tests for workflow.scripts.align_baker_scaffolds (cycle 03 sub-run A).

BAKER scaffolds carry the target as a fused chain B (HLA[1:180] + peptide);
the aligner superposes the HLA substring onto the truncated reference and
rewrites the complex into our A=binder / B=HLA / C=peptide layout.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio.PDB import PDBParser

from workflow.scripts import align_baker_scaffolds, align_scaffolds

REPO_ROOT = Path(__file__).resolve().parent.parent
MOCK_GLOB = "tests/fixtures/baker_truncated_mock/scaf*.pdb"
REFERENCE = REPO_ROOT / "tests/fixtures/targets/3hpj_baker_truncated_mock.pdb"


def _ca_counts(pdb: Path) -> dict[str, int]:
    model = PDBParser(QUIET=True).get_structure("x", str(pdb))[0]
    return {
        c.id: sum(1 for r in c.get_residues() if r.id[0].strip() == "" and "CA" in r)
        for c in model.get_chains()
    }


def test_aligned_layout_is_binder_hla_peptide(tmp_path: Path) -> None:
    aligned = align_baker_scaffolds.align_baker_scaffolds(
        MOCK_GLOB, REFERENCE, tmp_path / "aligned"
    )
    assert [p.name for p in aligned] == ["scaf0.pdb", "scaf1.pdb"]
    for p in aligned:
        counts = _ca_counts(p)
        assert set(counts) == {"A", "B", "C"}
        assert counts["B"] == 180  # reference HLA
        assert counts["C"] == 9  # reference peptide
        assert counts["A"] >= 3  # transformed binder


def test_hla_substring_superposition_is_recovered(tmp_path: Path) -> None:
    """Scaffold HLA is the reference HLA shifted +100 A in x, so superposing
    the first 180 CA collapses the RMSD to ~0."""
    rms = align_baker_scaffolds.align_one_baker(
        REPO_ROOT / "tests/fixtures/baker_truncated_mock/scaf0.pdb",
        REFERENCE,
        tmp_path / "scaf0.pdb",
    )
    assert rms < 1e-3


def test_missing_glob_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        align_baker_scaffolds.align_baker_scaffolds(
            "tests/fixtures/nope/scaf*.pdb", REFERENCE, tmp_path
        )


def test_main_mock_exits_zero(tmp_path: Path) -> None:
    rc = align_baker_scaffolds.main(["--mock", "--out-dir", str(tmp_path / "aligned")])
    assert rc == 0
    assert len(list((tmp_path / "aligned").glob("scaf*.pdb"))) == 2


def test_align_scaffolds_dispatches_to_baker_layout(tmp_path: Path) -> None:
    """align_scaffolds(..., baker_layout=True) must delegate and produce the
    A=binder / B=HLA / C=peptide layout, not the legacy passthrough."""
    aligned = align_scaffolds.align_scaffolds(
        MOCK_GLOB, REFERENCE, tmp_path / "aligned", baker_layout=True
    )
    assert len(aligned) == 2
    for p in aligned:
        assert set(_ca_counts(p)) == {"A", "B", "C"}
