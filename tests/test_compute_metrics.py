"""Unit tests for workflow.scripts.compute_metrics."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from workflow.scripts import compute_metrics


def _write_scores(path: Path, pae_key: str, pae_value: list[list[float]]) -> None:
    path.write_text(json.dumps({pae_key: pae_value, "plddt": [90.0]}))


def test_load_pae_new_key(tmp_path: Path) -> None:
    scores = tmp_path / "scores.json"
    _write_scores(scores, "pae", [[0.0, 1.0], [1.0, 0.0]])
    assert compute_metrics.load_pae(scores) == [[0.0, 1.0], [1.0, 0.0]]


def test_load_pae_legacy_key(tmp_path: Path) -> None:
    scores = tmp_path / "scores.json"
    _write_scores(scores, "predicted_aligned_error", [[2.0, 3.0], [3.0, 2.0]])
    assert compute_metrics.load_pae(scores) == [[2.0, 3.0], [3.0, 2.0]]


def test_load_pae_missing_raises(tmp_path: Path) -> None:
    scores = tmp_path / "scores.json"
    scores.write_text(json.dumps({"plddt": [90.0]}))
    with pytest.raises(KeyError, match="missing both"):
        compute_metrics.load_pae(scores)
