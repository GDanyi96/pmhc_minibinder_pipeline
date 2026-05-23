"""Unit tests for workflow.scripts.partial_diffuse (mock synthesis)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from workflow.scripts import partial_diffuse
from workflow.scripts.run_rfdiffusion import _geometry_pass, _hotspot_ca_xyz, _load_target

REPO_ROOT = Path(__file__).resolve().parent.parent
REFERENCE = REPO_ROOT / "tests/fixtures/stage1/mock_clean.pdb"
MANIFEST = REPO_ROOT / "tests/fixtures/stage1/mock_target.yaml"


def _write_ca_pdb(path: Path, chain_lengths: dict[str, int]) -> None:
    """Write a minimal CA-only PDB with the given per-chain residue counts."""
    lines: list[str] = []
    serial = 1
    for chain_id, n in chain_lengths.items():
        for resnum in range(1, n + 1):
            x, y, z = float(serial), 0.0, 0.0
            lines.append(
                f"ATOM  {serial:>5d}  CA  ALA {chain_id}{resnum:>4d}    "
                f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C"
            )
            serial += 1
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


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


def test_derive_contigs_subrun_a_80mer(tmp_path: Path) -> None:
    scaffold = tmp_path / "scaf0.pdb"
    _write_ca_pdb(scaffold, {"A": 80, "B": 180, "C": 9})
    assert partial_diffuse._derive_contigs_subrun_a(scaffold) == "[80-80/0 B1-180/0 C1-9]"


def test_derive_contigs_empty_chain_a_raises(tmp_path: Path) -> None:
    scaffold = tmp_path / "scaf_no_a.pdb"
    _write_ca_pdb(scaffold, {"B": 180, "C": 9})
    with pytest.raises(ValueError, match="not derivable"):
        partial_diffuse._derive_contigs_subrun_a(scaffold)


def test_partial_diffuse_one_passes_contigs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    scaffold = tmp_path / "scaf0.pdb"
    _write_ca_pdb(scaffold, {"A": 80, "B": 180, "C": 9})
    captured: dict[str, list[str]] = {}

    monkeypatch.setattr(partial_diffuse, "_resolve_rfdiffusion_python", lambda: "python")
    monkeypatch.setattr(
        partial_diffuse.subprocess,
        "run",
        lambda cmd, **kwargs: captured.update(cmd=cmd),  # type: ignore[arg-type]
    )

    partial_diffuse.partial_diffuse_one(
        input_pdb=scaffold,
        out_prefix=tmp_path / "design_3000",
        seed=3000,
        ckpt=tmp_path / "ckpt.pt",
        partial_T=15,
        noise_scale_ca=0,
        hydra_dir=tmp_path / "hydra",
    )
    assert "contigmap.contigs=[80-80/0 B1-180/0 C1-9]" in captured["cmd"]
