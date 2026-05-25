"""Unit tests for workflow.scripts.partial_diffuse (mock synthesis)."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from workflow.scripts import partial_diffuse
from workflow.scripts.run_rfdiffusion import _geometry_pass, _hotspot_ca_xyz, _load_target

REPO_ROOT = Path(__file__).resolve().parent.parent
REFERENCE = REPO_ROOT / "tests/fixtures/stage1/mock_clean.pdb"
MANIFEST = REPO_ROOT / "tests/fixtures/stage1/mock_target.yaml"


def test_synthesized_design_is_four_chain(tmp_path: Path) -> None:
    out = tmp_path / "design_99000.pdb"
    partial_diffuse._synthesize_mock_design(out, REFERENCE, binder_length=80)
    assert out.exists()
    from Bio.PDB import PDBParser

    model = PDBParser(QUIET=True).get_structure("d", str(out))[0]
    chains = {c.id for c in model.get_chains()}
    assert chains == {"A", "B", "C", "D"}
    binder = [r for r in model["D"].get_residues() if r.id[0].strip() == ""]
    assert len(binder) == 80


def test_synthesized_design_passes_geometry(tmp_path: Path) -> None:
    out = tmp_path / "design_99000.pdb"
    partial_diffuse._synthesize_mock_design(out, REFERENCE, binder_length=80)
    target = _load_target(MANIFEST).primary
    hotspot_xyz = _hotspot_ca_xyz(target, Path(target.cleaned_pdb))
    geom = _geometry_pass(out, hotspot_xyz, (70, 110), frozenset({"A", "B", "C"}))
    assert geom["geometry_pass"] is True
    assert geom["ca_contact_to_hotspots_n"] >= 3


def test_main_mock_exits_zero(tmp_path: Path) -> None:
    out = tmp_path / "design_99000.pdb"
    rc = partial_diffuse.main(
        ["--mock", "--input-pdb", str(REFERENCE), "--out-pdb", str(out), "--seed", "99000"]
    )
    assert rc == 0
    assert out.exists()


def test_default_partial_t_is_conservative() -> None:
    assert partial_diffuse.DEFAULT_PARTIAL_T == 10
    assert partial_diffuse.DEFAULT_NOISE_SCALE_CA == 0


def test_binder_backbone_chain_length() -> None:
    chain = partial_diffuse._binder_backbone_chain(12)
    coords = np.asarray([r["CA"].get_coord() for r in chain.get_residues()])
    assert coords.shape == (12, 3)
