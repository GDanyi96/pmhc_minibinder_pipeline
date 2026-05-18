"""Unit tests for workflow.scripts.run_rfdiffusion helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from workflow.scripts.prep_target import (
    BinderRange,
    ResolvedChain,
    ResolvedHotspot,
    ResolvedTarget,
)
from workflow.scripts.run_rfdiffusion import (
    SeedsConfig,
    _build_contigmap,
    _compute_seed,
    _geometry_pass,
)


def _write_pdb(path: Path, atoms: list[tuple[str, int, tuple[float, float, float]]]) -> None:
    """Write a minimal PDB with given (chain_id, resnum, xyz) Ca atoms."""
    lines: list[str] = []
    for i, (chain, resnum, xyz) in enumerate(atoms, start=1):
        x, y, z = xyz
        lines.append(
            f"ATOM  {i:>5d}  CA  ALA {chain}{resnum:>4d}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def _make_target(length_min: int = 70, length_max: int = 110) -> ResolvedTarget:
    return ResolvedTarget(
        target_id="wt1_a0201",
        pdb_id="3HPJ",
        cleaned_pdb="ignored",
        peptide_sequence="RMFPNAPYL",
        chains={
            "A": ResolvedChain(role="heavy_chain", length=275),
            "B": ResolvedChain(role="beta2m", length=99),
            "C": ResolvedChain(role="peptide", length=9),
            "D": ResolvedChain(role="binder"),
        },
        hotspots={
            "peptide": [
                ResolvedHotspot(token="C4", chain="C", resnum_clean=4, resnum_rcsb=4),
            ],
            "mhc": [
                ResolvedHotspot(token="A65", chain="A", resnum_clean=65, resnum_rcsb=65),
            ],
        },
        binder=BinderRange(length_min=length_min, length_max=length_max),
    )


def test_build_contigmap_reads_actual_lengths(tmp_path: Path) -> None:
    pdb = tmp_path / "mini_clean.pdb"
    atoms: list[tuple[str, int, tuple[float, float, float]]] = []
    for r in range(1, 6):  # chain A: 5 residues
        atoms.append(("A", r, (float(r), 0.0, 0.0)))
    for r in range(1, 4):  # chain B: 3 residues
        atoms.append(("B", r, (0.0, float(r), 0.0)))
    for r in range(1, 10):  # chain C: 9 residues
        atoms.append(("C", r, (0.0, 0.0, float(r))))
    _write_pdb(pdb, atoms)

    target = _make_target(length_min=70, length_max=110)
    contig = _build_contigmap(target, pdb)
    assert contig == "A1-5/0 B1-3/0 C1-9/0 70-110"


def test_build_contigmap_rejects_missing_chain(tmp_path: Path) -> None:
    pdb = tmp_path / "incomplete.pdb"
    _write_pdb(pdb, [("A", 1, (0.0, 0.0, 0.0))])  # only chain A
    target = _make_target()
    with pytest.raises(ValueError, match="chain B absent"):
        _build_contigmap(target, pdb)


SEEDS = SeedsConfig(
    formulas={"rfdiffusion": "cycle * 1000 + design_index"},
    reserved={
        "cycle_01_rfdiffusion": [1000, 1999],
        "cycle_02_rfdiffusion": [2000, 2199],
    },
)


def test_compute_seed_formula_cycle2_start() -> None:
    assert _compute_seed(2, 0, SEEDS) == 2000


def test_compute_seed_formula_cycle2_end() -> None:
    assert _compute_seed(2, 199, SEEDS) == 2199


def test_compute_seed_rejects_out_of_range() -> None:
    with pytest.raises(ValueError, match="outside reserved range"):
        _compute_seed(2, 200, SEEDS)


def test_compute_seed_rejects_unknown_cycle() -> None:
    with pytest.raises(ValueError, match="no reserved range"):
        _compute_seed(99, 0, SEEDS)


def test_compute_seed_rejects_unknown_formula() -> None:
    bad = SeedsConfig(
        formulas={"rfdiffusion": "cycle + design_index"},
        reserved={"cycle_02_rfdiffusion": [2000, 2199]},
    )
    with pytest.raises(ValueError, match="unsupported rfdiffusion seed formula"):
        _compute_seed(2, 0, bad)


HOTSPOTS_ON_AXIS = np.asarray(
    [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0]], dtype=np.float64
)


def _passing_binder_pdb(path: Path, n: int = 80) -> None:
    """80 Ca atoms placed along y axis starting near the hotspots."""
    atoms: list[tuple[str, int, tuple[float, float, float]]] = []
    for i in range(n):
        atoms.append(("A", i + 1, (5.0, float(i) * 3.9, 0.0)))
    _write_pdb(path, atoms)


def test_geometry_pass_three_contacts(tmp_path: Path) -> None:
    pdb = tmp_path / "design_pass.pdb"
    _passing_binder_pdb(pdb, n=80)
    result = _geometry_pass(pdb, HOTSPOTS_ON_AXIS, (70, 110))
    assert result["geometry_pass"] is True
    assert result["ca_contact_to_hotspots_n"] >= 3
    assert result["binder_length"] == 80
    assert result["internal_clash_count"] == 0


def test_geometry_pass_no_contacts(tmp_path: Path) -> None:
    pdb = tmp_path / "design_fail_distant.pdb"
    atoms = [("A", i + 1, (100.0, float(i) * 3.9, 0.0)) for i in range(80)]
    _write_pdb(pdb, atoms)
    result = _geometry_pass(pdb, HOTSPOTS_ON_AXIS, (70, 110))
    assert result["geometry_pass"] is False
    assert result["ca_contact_to_hotspots_n"] == 0
    assert "insufficient_hotspot_contacts" in result["fail_reasons"]


def test_geometry_pass_length_oob(tmp_path: Path) -> None:
    pdb = tmp_path / "design_fail_short.pdb"
    _passing_binder_pdb(pdb, n=50)
    result = _geometry_pass(pdb, HOTSPOTS_ON_AXIS, (70, 110))
    assert result["geometry_pass"] is False
    assert "length_out_of_range" in result["fail_reasons"]


def test_geometry_pass_internal_clash(tmp_path: Path) -> None:
    pdb = tmp_path / "design_fail_clash.pdb"
    atoms: list[tuple[str, int, tuple[float, float, float]]] = []
    for i in range(80):
        atoms.append(("A", i + 1, (5.0, float(i) * 3.9, 0.0)))
    atoms[40] = ("A", 41, (5.0, 0.0, 0.0))  # collides with atom 0
    _write_pdb(pdb, atoms)
    result = _geometry_pass(pdb, HOTSPOTS_ON_AXIS, (70, 110))
    assert result["geometry_pass"] is False
    assert result["internal_clash_count"] > 0
    assert "internal_ca_clash" in result["fail_reasons"]


def test_geometry_pass_empty_binder(tmp_path: Path) -> None:
    pdb = tmp_path / "design_empty.pdb"
    pdb.write_text("END\n")
    result = _geometry_pass(pdb, HOTSPOTS_ON_AXIS, (70, 110))
    assert result["geometry_pass"] is False
    assert result["binder_length"] == 0
