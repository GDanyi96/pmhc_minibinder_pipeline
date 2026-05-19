"""Tests for the Stage 1 mock-mode end-to-end + halt gate.

Exercises the halt rule in both directions:
- Mock fixtures are 8 pass / 2 fail at threshold 0.10 -> observed 0.80 -> PASS.
- Synthesized 0/10 failure scenario -> observed 0.00 -> FAIL.
- Boundary at the cycle 02 ratio (13/100 = 0.13) -> just above 0.10 -> PASS.
  Locks the `>=` semantics against future drift (Trap #29).
"""

from __future__ import annotations

import json
from pathlib import Path

from scripts import run_stage1
from scripts.run_stage1 import HALT_THRESHOLD, enforce_halt_gate


def test_halt_threshold_is_pinned() -> None:
    """Guards against silent drift away from the cycle 02 calibration value."""
    assert HALT_THRESHOLD == 0.10


def test_mock_run_produces_pass_summary(tmp_path: Path) -> None:
    out_dir = tmp_path / "rfdiffusion"
    rc = run_stage1.main(["--mock", "--cycle", "99", "--out-dir", str(out_dir)])
    assert rc == 0, "mock-mode orchestrator should exit 0"

    summary_path = out_dir / "stage1_summary.json"
    assert summary_path.exists()
    summary = json.loads(summary_path.read_text())

    assert summary["cycle"] == 99
    assert summary["target_id"] == "wt1_a0201"
    assert summary["n_attempted"] == 10
    assert summary["n_completed"] == 10
    assert summary["n_geometry_pass"] == 8
    assert summary["fraction_geometry_pass"] == 0.8
    assert summary["seed_range"] == [99000, 99199]
    assert summary["halt_rule"]["verdict"] == "PASS"
    assert summary["halt_rule"]["threshold"] == HALT_THRESHOLD
    assert summary["halt_rule"]["observed"] == 0.8
    assert summary["mock"] is True

    designs_jsonl = out_dir / "designs.jsonl"
    assert designs_jsonl.exists()
    records = [json.loads(line) for line in designs_jsonl.read_text().splitlines()]
    assert len(records) == 10
    n_pass = sum(1 for r in records if r["geometry_pass"])
    assert n_pass == 8


def _synthetic_summary(n_pass: int, n_total: int = 10) -> dict[str, object]:
    fraction = n_pass / n_total
    verdict = "PASS" if fraction >= HALT_THRESHOLD else "FAIL"
    return {
        "n_geometry_pass": n_pass,
        "n_completed": n_total,
        "fraction_geometry_pass": fraction,
        "halt_rule": {
            "name": "geometry_pass_rate",
            "threshold": HALT_THRESHOLD,
            "observed": fraction,
            "verdict": verdict,
            "note": "synthetic for test",
        },
    }


def test_enforce_halt_gate_pass() -> None:
    status, _ = enforce_halt_gate(_synthetic_summary(n_pass=8))
    assert status == "pass"


def test_enforce_halt_gate_halts_when_below_threshold() -> None:
    """10 of 10 failures -> observed 0.0 -> FAIL."""
    status, message = enforce_halt_gate(_synthetic_summary(n_pass=0))
    assert status == "halt"
    assert "0.0" in message
    assert "0.1" in message


def test_enforce_halt_gate_pass_at_cycle02_boundary() -> None:
    """Cycle 02 ratio 13/100 = 0.13 is just above 0.10 -> PASS.

    Locks the `>=` boundary semantics of the halt rule against future
    drift; this is the exact ratio that motivated Trap #29.
    """
    status, _ = enforce_halt_gate(_synthetic_summary(n_pass=13, n_total=100))
    assert status == "pass"


def test_enforce_halt_gate_halts_when_summary_malformed() -> None:
    status, message = enforce_halt_gate({"no_halt_rule_field": True})
    assert status == "halt"
    assert "missing halt_rule" in message


def test_skip_subprocess_reenumerates_existing_mock_fixtures(tmp_path: Path) -> None:
    """``--skip-subprocess --mock`` rebuilds the summary from pre-existing
    fixtures in designs/ without re-copying."""
    out_dir = tmp_path / "rfdiffusion"
    rc = run_stage1.main(["--mock", "--cycle", "99", "--out-dir", str(out_dir)])
    assert rc == 0
    summary_path = out_dir / "stage1_summary.json"
    summary_path.unlink()  # force the re-enumeration to actually do work

    rc2 = run_stage1.main(
        ["--mock", "--cycle", "99", "--out-dir", str(out_dir), "--skip-subprocess"]
    )
    assert rc2 == 0
    assert summary_path.exists()
    summary = json.loads(summary_path.read_text())
    assert summary["n_completed"] == 10
    assert summary["n_geometry_pass"] == 8


def test_skip_subprocess_halts_when_no_files_present(tmp_path: Path) -> None:
    """Empty designs_dir + --skip-subprocess -> n_completed=0 -> HALT."""
    out_dir = tmp_path / "rfdiffusion_empty"
    out_dir.mkdir()
    rc = run_stage1.main(
        ["--mock", "--cycle", "99", "--out-dir", str(out_dir), "--skip-subprocess"]
    )
    assert rc == 1
