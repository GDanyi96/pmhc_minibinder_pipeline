"""Tests for the cycle-03 Stage 1 sub-run + merge orchestrators (mock)."""

from __future__ import annotations

import json
from pathlib import Path

from scripts import merge_stage1_subruns, run_stage1_subrun


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
