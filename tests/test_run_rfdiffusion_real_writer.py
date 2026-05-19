"""Regression test for the writer/reader filename contract.

Stage 1's existing mock canary copies fixtures via ``_copy_mock_fixtures``,
which strips ``mock_`` from ``mock_design_NNNNN.pdb`` and produces files
that match the (then-broken) zero-padded enumerator expectation by
construction. The mock writer and reader were always trivially aligned.

Real RFdiffusion writes ``design_{seed}.pdb`` where seed comes from
``inference.design_startnum`` (= cycle*1000 + design_index for the
project's seed formula), yielding variable-digit names like
``design_2000.pdb``. The cycle 02 incident: 200 PDBs successfully
written, summary reported n_completed=0, pipeline HALTed.

This test synthesizes the real-writer file layout (seed-named PDBs in
designs/) and runs ``run_rfdiffusion(skip_subprocess=True, mock=False)``,
which exercises the reader against the real producer's naming
convention. It would have caught the cycle 02 drift.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import yaml

from workflow.scripts.run_rfdiffusion import run_rfdiffusion

CYCLE = 2
NUM_DESIGNS = 40
N_PASS = 4
FIXTURE_DIR = Path("tests/fixtures/stage1")


def _atom_line(serial: int, chain: str, resnum: int, xyz: tuple[float, float, float]) -> str:
    x, y, z = xyz
    return (
        f"ATOM  {serial:>5d}  CA  ALA {chain}{resnum:>4d}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C"
    )


def _write_pdb(path: Path, atoms: list[tuple[str, int, tuple[float, float, float]]]) -> None:
    lines = ["HEADER    SYNTHETIC REAL-WRITER FIXTURE"]
    for i, (chain, resnum, xyz) in enumerate(atoms, start=1):
        lines.append(_atom_line(i, chain, resnum, xyz))
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def _passing_backbone(
    n: int, chain: str = "D"
) -> list[tuple[str, int, tuple[float, float, float]]]:
    """Straight backbone through the mock hotspot centroid (matches
    ``_make_fixtures._passing_backbone``).

    Default chain ``D`` matches mock_target.yaml's binder declaration so
    that ``_binder_ca_coords`` (which now identifies the binder by
    exclusion from motif chains A/B/C) selects this chain correctly.
    """
    atoms: list[tuple[str, int, tuple[float, float, float]]] = []
    for i in range(n):
        x = 1.5
        y = i * 3.9 - (n - 1) * 3.9 / 2 + 3.5
        z = 2.5
        atoms.append((chain, i + 1, (x, y, z)))
    return atoms


def _distant_backbone(
    n: int, chain: str = "D"
) -> list[tuple[str, int, tuple[float, float, float]]]:
    """Geometry-fail: 100 A from any hotspot.

    Default chain ``D`` matches mock_target.yaml's binder declaration.
    """
    atoms: list[tuple[str, int, tuple[float, float, float]]] = []
    for i in range(n):
        atoms.append((chain, i + 1, (100.0, float(i) * 3.9, 100.0)))
    return atoms


def _straight_motif_stub(
    chain: str, n: int, x_offset: float
) -> list[tuple[str, int, tuple[float, float, float]]]:
    """Filler atoms for a motif chain. Placed far from the hotspot cluster
    so they do not contribute Ca-contact hits; only their chain identity
    matters for ``_binder_ca_coords`` exclusion logic.
    """
    atoms: list[tuple[str, int, tuple[float, float, float]]] = []
    for i in range(n):
        atoms.append((chain, i + 1, (x_offset, float(i) * 3.9, 50.0)))
    return atoms


def _seeds_yaml(tmp_path: Path) -> Path:
    """seeds.yaml that reserves cycle 02's seed range [2000, 2999]."""
    path = tmp_path / "seeds.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "formulas": {"rfdiffusion": "cycle * 1000 + design_index"},
                "reserved": {"cycle_02_rfdiffusion": [2000, 2999]},
            }
        )
    )
    return path


def _rfdiff_yaml(tmp_path: Path) -> Path:
    """Real-mode rfdiffusion config; ckpt_override_path needs to exist
    on disk for ``_file_sha256`` even though we never invoke the subprocess."""
    ckpt = tmp_path / "fake_ckpt.pt"
    ckpt.write_bytes(b"")
    path = tmp_path / "rfdiff.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "inference": {
                    "num_designs": NUM_DESIGNS,
                    "ckpt_override_path": str(ckpt),
                },
            }
        )
    )
    return path


def _seed_for(design_index: int) -> int:
    return CYCLE * 1000 + design_index


def test_real_writer_enumeration_finds_seed_named_files(tmp_path: Path) -> None:
    """Synthesize the real-writer file layout, then run the enumerator.

    Pre-cycle-02-fix this fails with ``n_completed == 0`` because the
    reader expected ``design_{design_index:05d}.pdb`` while the
    synthesizer wrote ``design_{seed}.pdb`` (seed = cycle*1000 +
    design_index, matching RFdiffusion's inference.design_startnum).
    """
    out_dir = tmp_path / "rfdiffusion"
    designs_dir = out_dir / "designs"
    designs_dir.mkdir(parents=True)

    for design_index in range(NUM_DESIGNS):
        seed = _seed_for(design_index)
        path = designs_dir / f"design_{seed}.pdb"
        if design_index < N_PASS:
            _write_pdb(path, _passing_backbone(85))
        else:
            _write_pdb(path, _distant_backbone(85))

    target_manifest = FIXTURE_DIR / "mock_target.yaml"

    summary_path = run_rfdiffusion(
        target_manifest=target_manifest,
        rfdiff_yaml=_rfdiff_yaml(tmp_path),
        seeds_yaml=_seeds_yaml(tmp_path),
        cycle=CYCLE,
        out_dir=out_dir,
        mock=False,
        skip_subprocess=True,
    )
    summary: dict[str, Any] = json.loads(summary_path.read_text())

    assert summary["n_attempted"] == NUM_DESIGNS
    assert summary["n_completed"] == NUM_DESIGNS
    assert summary["n_geometry_pass"] == N_PASS
    assert summary["fraction_geometry_pass"] == N_PASS / NUM_DESIGNS
    assert summary["mock"] is False
    assert summary["seed_range"] == [2000, 2999]

    designs_jsonl = out_dir / "designs.jsonl"
    records = [json.loads(line) for line in designs_jsonl.read_text().splitlines()]
    assert len(records) == NUM_DESIGNS
    assert all(r["design_id"].startswith("design_2") for r in records)
    assert not (designs_dir / "design_00000.pdb").exists()


def test_real_writer_missing_files_counted_as_incomplete(tmp_path: Path) -> None:
    """If RFdiffusion crashed mid-run, only the seeds that landed on disk
    should count toward n_completed; n_attempted still equals num_designs."""
    out_dir = tmp_path / "rfdiffusion"
    designs_dir = out_dir / "designs"
    designs_dir.mkdir(parents=True)

    n_present = 12
    for design_index in range(n_present):
        seed = _seed_for(design_index)
        _write_pdb(designs_dir / f"design_{seed}.pdb", _passing_backbone(85))

    summary_path = run_rfdiffusion(
        target_manifest=FIXTURE_DIR / "mock_target.yaml",
        rfdiff_yaml=_rfdiff_yaml(tmp_path),
        seeds_yaml=_seeds_yaml(tmp_path),
        cycle=CYCLE,
        out_dir=out_dir,
        mock=False,
        skip_subprocess=True,
    )
    summary = json.loads(summary_path.read_text())

    assert summary["n_attempted"] == NUM_DESIGNS
    assert summary["n_completed"] == n_present


def test_real_writer_binder_chain_is_D(tmp_path: Path) -> None:
    """Synthesize the real RFdiffusion chain layout: motif on A/B/C,
    binder on D. Pre-fix ``_binder_ca_coords`` selected chain A (the
    275-res HLA heavy chain) and tripped length-out-of-range on 200/200
    designs in cycle 02. After the fix, the binder is identified by
    exclusion (chains - fixed_chains) and reports the correct length.

    Reference: docs/known_traps.md trap #28.
    """
    out_dir = tmp_path / "rfdiffusion"
    designs_dir = out_dir / "designs"
    designs_dir.mkdir(parents=True)

    n_designs_here = 4
    for design_index in range(n_designs_here):
        seed = _seed_for(design_index)
        path = designs_dir / f"design_{seed}.pdb"
        atoms: list[tuple[str, int, tuple[float, float, float]]] = []
        atoms.extend(_straight_motif_stub("A", 275, x_offset=80.0))
        atoms.extend(_straight_motif_stub("B", 100, x_offset=120.0))
        atoms.extend(_straight_motif_stub("C", 9, x_offset=160.0))
        atoms.extend(_passing_backbone(85, chain="D"))
        _write_pdb(path, atoms)

    rfdiff_yaml_path = _rfdiff_yaml(tmp_path)
    # Override num_designs to match what we wrote.
    cfg = yaml.safe_load(rfdiff_yaml_path.read_text())
    cfg["inference"]["num_designs"] = n_designs_here
    rfdiff_yaml_path.write_text(yaml.safe_dump(cfg))

    summary_path = run_rfdiffusion(
        target_manifest=FIXTURE_DIR / "mock_target.yaml",
        rfdiff_yaml=rfdiff_yaml_path,
        seeds_yaml=_seeds_yaml(tmp_path),
        cycle=CYCLE,
        out_dir=out_dir,
        mock=False,
        skip_subprocess=True,
    )
    summary: dict[str, Any] = json.loads(summary_path.read_text())

    designs_jsonl = out_dir / "designs.jsonl"
    records = [json.loads(line) for line in designs_jsonl.read_text().splitlines()]
    assert len(records) == n_designs_here
    for record in records:
        assert record["binder_length"] == 85, record
    assert summary["n_geometry_pass"] == n_designs_here


def test_binder_ca_coords_raises_on_ambiguous_chains(tmp_path: Path) -> None:
    """When more than one candidate chain remains after subtracting the
    declared motif chains, ``_binder_ca_coords`` must raise ValueError.
    This is the loud-failure guarantee for the chain-exclusion contract.

    Reference: docs/known_traps.md trap #28.
    """
    from workflow.scripts.run_rfdiffusion import _binder_ca_coords

    pdb_path = tmp_path / "ambiguous.pdb"
    atoms: list[tuple[str, int, tuple[float, float, float]]] = []
    atoms.extend(_straight_motif_stub("A", 10, x_offset=0.0))
    atoms.extend(_straight_motif_stub("B", 10, x_offset=20.0))
    atoms.extend(_straight_motif_stub("C", 10, x_offset=40.0))
    atoms.extend(_passing_backbone(20, chain="D"))
    _write_pdb(pdb_path, atoms)

    import pytest

    with pytest.raises(ValueError) as exc:
        _binder_ca_coords(pdb_path, fixed_chains=frozenset({"A"}))

    message = str(exc.value)
    assert "candidates" in message
    assert "['A', 'B', 'C', 'D']" in message or "'B'" in message
