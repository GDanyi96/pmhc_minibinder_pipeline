"""Tests for the cycle-03 Stage 1 sub-run + merge orchestrators (mock)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts import merge_stage1_subruns, run_stage1_subrun
from workflow.scripts.run_rfdiffusion import _load_target

REPO_ROOT = Path(__file__).resolve().parent.parent
MOCK_MANIFEST = REPO_ROOT / "tests/fixtures/stage1/mock_target.yaml"
MOCK_REFERENCE = REPO_ROOT / "tests/fixtures/stage1/mock_clean.pdb"


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


def test_fixed_chains_subrun_a_real_is_truncated_motif() -> None:
    target = _load_target(MOCK_MANIFEST).primary
    assert run_stage1_subrun._fixed_chains_for_subrun("a", target, mock=False) == frozenset(
        {"B", "C"}
    )
    manifest_set = frozenset(c for c, ch in target.chains.items() if ch.role != "binder")
    assert run_stage1_subrun._fixed_chains_for_subrun("b", target, mock=False) == manifest_set
    assert run_stage1_subrun._fixed_chains_for_subrun("a", target, mock=True) == manifest_set


def test_subrun_a_real_scores_three_chain_design(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """End-to-end call-site guard (Trap #28 lesson): drive the real run() path
    with GPU/IO boundaries stubbed and a 3-chain truncated design (A=binder,
    B=HLA, C=peptide). The binder must be identified through the real
    _fixed_chains_for_subrun -> _geometry_pass -> _binder_ca_coords chain, not
    crash on the empty-candidate exclusion."""
    aligned = tmp_path / "aligned" / "scaf0.pdb"
    aligned.parent.mkdir(parents=True)
    _write_ca_pdb(aligned, {"A": 80, "B": 180, "C": 9})

    def fake_align(scaffold_glob, reference_pdb, out_dir, **kwargs):  # type: ignore[no-untyped-def]
        return [aligned]

    def fake_diffuse(input_pdb, out_prefix, **kwargs):  # type: ignore[no-untyped-def]
        _write_ca_pdb(Path(out_prefix).with_suffix(".pdb"), {"A": 80, "B": 180, "C": 9})

    monkeypatch.setattr(run_stage1_subrun.align_scaffolds, "align_scaffolds", fake_align)
    monkeypatch.setattr(run_stage1_subrun.partial_diffuse, "partial_diffuse_one", fake_diffuse)
    monkeypatch.setattr(run_stage1_subrun, "_file_sha256", lambda p: "deadbeef")
    monkeypatch.setattr(run_stage1_subrun, "_git_sha", lambda p: "abc1234")

    rc = run_stage1_subrun.run(
        subrun="a",
        config_yaml=REPO_ROOT / "configs/rfdiffusion_subrun_a.yaml",
        seeds_yaml=REPO_ROOT / "configs/seeds.yaml",
        target_manifest=MOCK_MANIFEST,
        reference_pdb=MOCK_REFERENCE,
        cycle=99,
        out_dir=tmp_path / "subrun_a",
        mock=False,
    )
    assert rc == 0
    records = [
        json.loads(line)
        for line in (tmp_path / "subrun_a" / "designs.jsonl").read_text().splitlines()
        if line.strip()
    ]
    assert len(records) == 1
    assert records[0]["binder_length"] == 80


def test_subrun_a_writer_reader_filename_roundtrip(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Trap #33 contract test: drive the REAL partial_diffuse_one (prefix
    construction) and emulate RFdiffusion's filename rule (it appends
    ``_{design_startnum}`` to ``inference.output_prefix``), then let run()'s
    real enumeration glob find the design. Reintroducing the seed into the
    prefix would yield ``design_{seed}_{seed}.pdb`` -> n_completed: 0, the exact
    symptom that was invisible until the symlink recovery."""
    aligned = tmp_path / "aligned" / "scaf0.pdb"
    aligned.parent.mkdir(parents=True)
    _write_ca_pdb(aligned, {"A": 71, "B": 180, "C": 9})

    def fake_align(scaffold_glob, reference_pdb, out_dir, **kwargs):  # type: ignore[no-untyped-def]
        return [aligned]

    def fake_rfdiffusion_run(cmd, **kwargs):  # type: ignore[no-untyped-def]
        prefix = next(
            a.split("=", 1)[1] for a in cmd if a.startswith("inference.output_prefix=")
        )
        startnum = next(
            a.split("=", 1)[1] for a in cmd if a.startswith("inference.design_startnum=")
        )
        # RFdiffusion appends _{startnum} to the prefix.
        _write_ca_pdb(Path(f"{prefix}_{startnum}.pdb"), {"A": 71, "B": 180, "C": 9})

    monkeypatch.setattr(run_stage1_subrun.align_scaffolds, "align_scaffolds", fake_align)
    monkeypatch.setattr(run_stage1_subrun.partial_diffuse, "_resolve_rfdiffusion_python", lambda: "python")
    monkeypatch.setattr(run_stage1_subrun.partial_diffuse.subprocess, "run", fake_rfdiffusion_run)
    monkeypatch.setattr(run_stage1_subrun, "_file_sha256", lambda p: "deadbeef")
    monkeypatch.setattr(run_stage1_subrun, "_git_sha", lambda p: "abc1234")

    rc = run_stage1_subrun.run(
        subrun="a",
        config_yaml=REPO_ROOT / "configs/rfdiffusion_subrun_a.yaml",
        seeds_yaml=REPO_ROOT / "configs/seeds.yaml",
        target_manifest=MOCK_MANIFEST,
        reference_pdb=MOCK_REFERENCE,
        cycle=99,
        out_dir=tmp_path / "subrun_a",
        mock=False,
    )
    assert rc == 0
    summary = json.loads((tmp_path / "subrun_a" / "subrun_summary.json").read_text())
    assert summary["n_completed"] == 1
    # The single design lands as a single-suffix file the enumerator found.
    designs = list((tmp_path / "subrun_a" / "designs").glob("design_*.pdb"))
    assert [p.name for p in designs] == ["design_99000.pdb"]


def _run_subrun(subrun: str, out_dir: Path) -> dict:
    rc = run_stage1_subrun.main(
        ["--mock", "--subrun", subrun, "--cycle", "99", "--out-dir", str(out_dir)]
    )
    assert rc == 0
    return json.loads((out_dir / "subrun_summary.json").read_text())


def test_subrun_a_aligns_and_diffuses_each_scaffold(tmp_path: Path) -> None:
    summary = _run_subrun("a", tmp_path / "subrun_a")
    assert summary["subrun"] == "a"
    assert summary["n_completed"] == 3  # one per mock scaffold
    assert summary["n_geometry_pass"] == 3
    assert summary["partial_T"] == 15
    assert (tmp_path / "subrun_a" / "aligned").is_dir()


def test_subrun_b_diffuses_from_seed(tmp_path: Path) -> None:
    summary = _run_subrun("b", tmp_path / "subrun_b")
    assert summary["subrun"] == "b"
    assert summary["n_completed"] == 2
    assert summary["partial_T"] == 10
    assert summary["seed_offset"] == 150


def test_subrun_seeds_are_disjoint_and_in_range(tmp_path: Path) -> None:
    _run_subrun("a", tmp_path / "a")
    _run_subrun("b", tmp_path / "b")
    seeds_a = {
        json.loads(line)["seed"]
        for line in (tmp_path / "a" / "designs.jsonl").read_text().splitlines()
    }
    seeds_b = {
        json.loads(line)["seed"]
        for line in (tmp_path / "b" / "designs.jsonl").read_text().splitlines()
    }
    assert seeds_a.isdisjoint(seeds_b)
    assert all(99000 <= s <= 99199 for s in seeds_a | seeds_b)


def test_merge_combines_subruns_and_passes(tmp_path: Path) -> None:
    _run_subrun("a", tmp_path / "subrun_a")
    _run_subrun("b", tmp_path / "subrun_b")
    rc = merge_stage1_subruns.main(
        [
            "--mock",
            "--cycle",
            "99",
            "--subrun-dir",
            str(tmp_path / "subrun_a"),
            "--subrun-dir",
            str(tmp_path / "subrun_b"),
            "--out-dir",
            str(tmp_path / "rfdiffusion"),
        ]
    )
    assert rc == 0
    summary = json.loads((tmp_path / "rfdiffusion" / "stage1_summary.json").read_text())
    assert summary["n_completed"] == 5
    assert summary["halt_rule"]["verdict"] == "PASS"
    designs = list((tmp_path / "rfdiffusion" / "designs").glob("design_*.pdb"))
    assert len(designs) == 5
    # designs.jsonl is the canonical handoff Stage 2 reads.
    records = (tmp_path / "rfdiffusion" / "designs.jsonl").read_text().splitlines()
    assert len(records) == 5


def test_merge_halts_when_all_fail_geometry(tmp_path: Path) -> None:
    subrun_dir = tmp_path / "subrun_bad"
    (subrun_dir / "designs").mkdir(parents=True)
    with (subrun_dir / "designs.jsonl").open("w") as fh:
        for seed in (99000, 99001):
            fh.write(
                json.dumps(
                    {
                        "design_id": f"design_{seed}",
                        "pdb_path": str(subrun_dir / "designs" / f"design_{seed}.pdb"),
                        "target_id": "wt1_a0201",
                        "binder_length": 80,
                        "geometry_pass": False,
                    }
                )
                + "\n"
            )
    rc = merge_stage1_subruns.main(
        [
            "--mock",
            "--cycle",
            "99",
            "--subrun-dir",
            str(subrun_dir),
            "--out-dir",
            str(tmp_path / "rfdiffusion"),
        ]
    )
    assert rc == 1
    summary = json.loads((tmp_path / "rfdiffusion" / "stage1_summary.json").read_text())
    assert summary["halt_rule"]["verdict"] == "FAIL"
