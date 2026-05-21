"""Unit tests for workflow.scripts.contact_filter (cycle 03 placement gate)."""

from __future__ import annotations

import json
from pathlib import Path

from workflow.scripts import contact_filter, partial_diffuse

REPO_ROOT = Path(__file__).resolve().parent.parent
REFERENCE = REPO_ROOT / "tests/fixtures/stage1/mock_clean.pdb"
SEED = REPO_ROOT / "tests/fixtures/design_2079_mock_seed.pdb"


def _make_design(tmp_path: Path) -> Path:
    out = tmp_path / "design.pdb"
    partial_diffuse._synthesize_mock_design(out, REFERENCE, binder_length=80)
    return out


def test_centroid_binder_contacts_peptide(tmp_path: Path) -> None:
    design = _make_design(tmp_path)
    n = contact_filter.count_peptide_contacts(design)
    assert n >= contact_filter.DEFAULT_MIN_CONTACTS


def test_seed_complex_contacts_peptide() -> None:
    """The stitched seed (chain D through the cluster) contacts the peptide."""
    assert contact_filter.passes(SEED) is True


def test_high_threshold_fails(tmp_path: Path) -> None:
    design = _make_design(tmp_path)
    assert contact_filter.passes(design, min_contacts=10_000) is False


def test_distance_monotonic(tmp_path: Path) -> None:
    design = _make_design(tmp_path)
    near = contact_filter.count_peptide_contacts(design, distance=3.0)
    far = contact_filter.count_peptide_contacts(design, distance=12.0)
    assert far >= near


def test_main_mock_runs_on_seed_complex(tmp_path: Path) -> None:
    out = tmp_path / "contacts.json"
    rc = contact_filter.main(["--mock", "--out", str(out)])
    assert rc == 0
    records = json.loads(out.read_text())
    assert len(records) == 1
    assert records[0]["pass"] is True


def test_main_writes_records_for_design_set(tmp_path: Path) -> None:
    for i in range(3):
        partial_diffuse._synthesize_mock_design(
            tmp_path / f"design_{i:05d}.pdb", REFERENCE, binder_length=80
        )
    out = tmp_path / "contacts.json"
    rc = contact_filter.main(["--input-glob", str(tmp_path / "design_*.pdb"), "--out", str(out)])
    assert rc == 0
    records = json.loads(out.read_text())
    assert len(records) == 3
    assert all("peptide_contacts" in r and "pass" in r for r in records)
