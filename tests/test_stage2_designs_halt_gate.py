"""Tests for the Stage 2 designs mock-mode orchestrator + halt gate.

Mirrors tests/test_stage1_halt_gate.py. The fixture is calibrated so that
exactly 4 of 40 mock AF2 predictions pass the intermediate cut, i.e.
observed = 0.10 exactly at the >= boundary -- a refactor canary for the
halt rule's PASS-at-boundary semantics.
"""

from __future__ import annotations

import json
from pathlib import Path

from scripts import run_stage2
from scripts.run_stage2 import _HALT_THRESHOLD, enforce_halt_gate


def test_halt_threshold_is_pinned() -> None:
    """Guards against silent drift away from the calibration-only 0.10."""
    assert _HALT_THRESHOLD == 0.10


def test_mock_run_produces_pass_summary(tmp_path: Path) -> None:
    out_dir = tmp_path / "stage2"
    rc = run_stage2.main(["--mock", "--cycle", "99", "--out-dir", str(out_dir)])
    assert rc == 0, "mock-mode orchestrator should exit 0"

    summary_path = out_dir / "stage2_summary.json"
    assert summary_path.exists()
    summary = json.loads(summary_path.read_text())

    assert summary["cycle"] == 99
    assert summary["target_id"] == "wt1_a0201"
    assert summary["n_backbones_in"] == 10
    assert summary["n_sequences_designed"] == 40
    assert summary["n_af2_folded"] == 40
    assert summary["n_pass_intermediate"] == 4
    assert summary["fraction_pass_intermediate"] == 0.10
    assert summary["halt_rule"]["verdict"] == "PASS"
    assert summary["halt_rule"]["threshold"] == _HALT_THRESHOLD
    assert summary["halt_rule"]["observed"] == 0.10
    assert summary["mock"] is True

    sequences_jsonl = out_dir / "proteinmpnn_designs" / "sequences.jsonl"
    assert sequences_jsonl.exists()
    records = [json.loads(line) for line in sequences_jsonl.read_text().splitlines()]
    assert len(records) == 40

    metrics_jsonl = out_dir / "af2_designs" / "metrics.jsonl"
    assert metrics_jsonl.exists()
    metric_records = [json.loads(line) for line in metrics_jsonl.read_text().splitlines()]
    assert len(metric_records) == 40


def _synthetic_summary(n_pass: int, n_total: int = 40) -> dict[str, object]:
    fraction = n_pass / n_total
    verdict = "PASS" if fraction >= _HALT_THRESHOLD else "FAIL"
    return {
        "n_pass_intermediate": n_pass,
        "n_af2_folded": n_total,
        "fraction_pass_intermediate": fraction,
        "halt_rule": {
            "name": "fraction_pass_intermediate",
            "threshold": _HALT_THRESHOLD,
            "observed": fraction,
            "verdict": verdict,
            "note": "synthetic",
        },
    }


def test_enforce_halt_gate_pass_at_boundary() -> None:
    """4/40 = 0.10 exactly -> PASS (refactor canary for >= semantics)."""
    status, _ = enforce_halt_gate(_synthetic_summary(n_pass=4))
    assert status == "pass"


def test_enforce_halt_gate_halts_when_below_threshold() -> None:
    """3/40 = 0.075 < 0.10 -> FAIL."""
    status, message = enforce_halt_gate(_synthetic_summary(n_pass=3))
    assert status == "halt"
    assert "0.075" in message
    assert "0.1" in message


def test_enforce_halt_gate_halts_when_summary_malformed() -> None:
    status, message = enforce_halt_gate({"no_halt_rule_field": True})
    assert status == "halt"
    assert "missing halt_rule" in message


def test_metrics_only_replays_existing_summary(tmp_path: Path) -> None:
    out_dir = tmp_path / "stage2"
    rc = run_stage2.main(["--mock", "--cycle", "99", "--out-dir", str(out_dir)])
    assert rc == 0
    rc2 = run_stage2.main(["--mock", "--cycle", "99", "--out-dir", str(out_dir), "--metrics-only"])
    assert rc2 == 0


def test_metrics_only_errors_when_summary_missing(tmp_path: Path) -> None:
    out_dir = tmp_path / "stage2_empty"
    out_dir.mkdir()
    rc = run_stage2.main(["--mock", "--cycle", "99", "--out-dir", str(out_dir), "--metrics-only"])
    assert rc == 1


def test_fan_in_top_n_mock_caps_at_record_count(tmp_path: Path) -> None:
    """In mock mode the orchestrator caps fan-in at min(40, n_records)
    rather than the real-mode value of 100 -- ensures the boundary fixture
    is exercised end-to-end without scaling up the AF2 fixture set."""
    out_dir = tmp_path / "stage2"
    rc = run_stage2.main(["--mock", "--cycle", "99", "--out-dir", str(out_dir)])
    assert rc == 0
    summary = json.loads((out_dir / "stage2_summary.json").read_text())
    assert summary["n_af2_folded"] == 40
