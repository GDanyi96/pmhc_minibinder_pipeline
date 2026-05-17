"""Tests for scripts.run_controls.enforce_halt_gate (commit-5 rule rewrite)."""

from __future__ import annotations

import pytest

from scripts import run_controls
from scripts.run_controls import ControlsThresholds, enforce_halt_gate

DEFAULT = ControlsThresholds(
    p1_ipae_max=7.0,
    p2_ipae_tolerance=2.0,
    pos_neg_ipae_gap_min=10.0,
    pos_neg_iplddt_gap_min=30.0,
)


def _metrics(ipae: dict[str, float], iplddt: dict[str, float]) -> dict:
    return {cid: {"ipae": ipae[cid], "iplddt": iplddt[cid]} for cid in ipae}


POD_CYCLE1_IPAE = {"P1": 4.88, "P2": 4.56, "P3": 4.54, "N1": 24.72, "N2": 25.70}
POD_CYCLE1_IPLDDT = {"P1": 94.50, "P2": 95.95, "P3": 94.21, "N1": 34.48, "N2": 31.43}


def test_pod_cycle1_passes_new_halt_gate() -> None:
    """The literal cycle 1 pod numbers must pass the new rule.

    The previous P1-rank rule fired false on these. The new rule (positives
    vs negatives gap) is the de-risking event for Stage 1.
    """
    status, msg = enforce_halt_gate(_metrics(POD_CYCLE1_IPAE, POD_CYCLE1_IPLDDT), DEFAULT)
    assert status == "pass", msg


def test_halt_when_p1_above_threshold() -> None:
    ipae = {"P1": 9.0, "P2": 6.5, "P3": 5.0, "N1": 20.0, "N2": 22.0}
    iplddt = {"P1": 90.0, "P2": 92.0, "P3": 91.0, "N1": 30.0, "N2": 30.0}
    status, msg = enforce_halt_gate(_metrics(ipae, iplddt), DEFAULT)
    assert status == "halt"
    assert "P1 iPAE" in msg


def test_halt_when_p2_calibration_off() -> None:
    ipae = {"P1": 5.0, "P2": 10.0, "P3": 5.0, "N1": 25.0, "N2": 27.0}
    iplddt = {"P1": 90.0, "P2": 92.0, "P3": 91.0, "N1": 30.0, "N2": 30.0}
    status, msg = enforce_halt_gate(_metrics(ipae, iplddt), DEFAULT)
    assert status == "halt"
    assert "P2" in msg and "calibration" in msg


def test_halt_when_pos_neg_ipae_gap_too_small() -> None:
    ipae = {"P1": 5.0, "P2": 6.5, "P3": 5.0, "N1": 12.0, "N2": 14.0}
    iplddt = {"P1": 90.0, "P2": 92.0, "P3": 91.0, "N1": 30.0, "N2": 30.0}
    status, msg = enforce_halt_gate(_metrics(ipae, iplddt), DEFAULT)
    assert status == "halt"
    assert "iPAE gap" in msg


def test_halt_when_pos_neg_iplddt_gap_too_small() -> None:
    ipae = {"P1": 5.0, "P2": 6.5, "P3": 5.0, "N1": 25.0, "N2": 27.0}
    iplddt = {"P1": 70.0, "P2": 72.0, "P3": 71.0, "N1": 50.0, "N2": 50.0}
    status, msg = enforce_halt_gate(_metrics(ipae, iplddt), DEFAULT)
    assert status == "halt"
    assert "ipLDDT gap" in msg


def test_specificity_note_logged_at_info(caplog: pytest.LogCaptureFixture) -> None:
    """The P3 specificity note must log at INFO, not WARNING."""
    import logging

    with caplog.at_level(logging.INFO, logger=run_controls.logger.name):
        enforce_halt_gate(_metrics(POD_CYCLE1_IPAE, POD_CYCLE1_IPLDDT), DEFAULT)
    info_records = [r for r in caplog.records if r.levelno == logging.INFO]
    assert any("specificity" in r.message for r in info_records)
    assert not any(r.levelno == logging.WARNING for r in caplog.records)
