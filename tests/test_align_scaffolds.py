"""Unit tests for workflow.scripts.align_scaffolds (cycle 03 sub-run A)."""

from __future__ import annotations

from pathlib import Path

import pytest

from workflow.scripts import align_scaffolds

REPO_ROOT = Path(__file__).resolve().parent.parent
MOCK_GLOB = "tests/fixtures/baker_library_mock/scaf*.pdb"
REFERENCE = REPO_ROOT / "tests/fixtures/stage1/mock_clean.pdb"


def test_aligns_all_mock_scaffolds(tmp_path: Path) -> None:
    aligned = align_scaffolds.align_scaffolds(MOCK_GLOB, REFERENCE, tmp_path / "aligned")
    assert len(aligned) == 3
    assert all(p.exists() for p in aligned)
    assert [p.name for p in aligned] == ["scaf0.pdb", "scaf1.pdb", "scaf2.pdb"]


def test_rigid_translation_is_recovered(tmp_path: Path) -> None:
    """Scaffold target chain B is the reference MHC chain A shifted +100 A in
    x, so the superposition RMSD must collapse to ~0."""
    rms = align_scaffolds.align_one(
        REPO_ROOT / "tests/fixtures/baker_library_mock/scaf0.pdb",
        REFERENCE,
        tmp_path / "scaf0.pdb",
    )
    assert rms < 1e-3


def test_missing_glob_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        align_scaffolds.align_scaffolds("tests/fixtures/nope/scaf*.pdb", REFERENCE, tmp_path)


def test_main_mock_exits_zero(tmp_path: Path) -> None:
    rc = align_scaffolds.main(["--mock", "--out-dir", str(tmp_path / "aligned")])
    assert rc == 0
    assert len(list((tmp_path / "aligned").glob("scaf*.pdb"))) == 3
