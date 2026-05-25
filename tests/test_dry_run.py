"""Smoke test: the pipeline DAG resolves in mock mode."""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not installed")
def test_snakemake_dry_run_mock() -> None:
    """`snakemake --dry-run --config mock=true -j1` must exit 0."""
    result = subprocess.run(
        ["snakemake", "--dry-run", "--config", "mock=true", "-j1"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert (
        result.returncode == 0
    ), f"snakemake dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not installed")
def test_snakemake_mock_end_to_end(tmp_path: Path) -> None:
    """`snakemake --config mock=true -j1` runs the full DAG via fixture cp."""
    work = tmp_path / "run"
    shutil.copytree(REPO_ROOT, work, ignore=shutil.ignore_patterns(".git", "results"))
    result = subprocess.run(
        ["snakemake", "--config", "mock=true", "-j1", "--quiet"],
        cwd=work,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert (
        result.returncode == 0
    ), f"snakemake mock run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    assert (work / "reports" / "cycle_01.md").exists()


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not installed")
def test_snakemake_dry_run_cycle03_partial() -> None:
    """Cycle 03 (>=3) resolves the partial-diffusion DAG (subruns + merge)."""
    result = subprocess.run(
        ["snakemake", "--dry-run", "--config", "mock=true", "cycle=03", "-j1"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, f"cycle03 dry-run failed:\n{result.stdout}\n{result.stderr}"
    assert "stage1_merge" in result.stdout
    assert "stage1_subrun_a" in result.stdout
    assert "stage1_subrun_b" in result.stdout


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not installed")
def test_snakemake_mock_end_to_end_cycle99(tmp_path: Path) -> None:
    """Full mock pipeline at cycle 99 exercises the partial-diffusion path
    (sub-runs A+B -> merge -> Stage 2 -> report) end to end, no GPU."""
    work = tmp_path / "run"
    shutil.copytree(REPO_ROOT, work, ignore=shutil.ignore_patterns(".git", "results", ".snakemake"))
    result = subprocess.run(
        ["snakemake", "--config", "mock=true", "cycle=99", "-j1", "--quiet"],
        cwd=work,
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert (
        result.returncode == 0
    ), f"cycle99 mock run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    assert (work / "reports" / "cycle_99.md").exists()
    merged = work / "results" / "cycle_99" / "stage1" / "rfdiffusion" / "stage1_summary.json"
    assert json.loads(merged.read_text())["halt_rule"]["verdict"] == "PASS"
