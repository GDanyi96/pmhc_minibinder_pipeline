"""Smoke test: the pipeline DAG resolves in mock mode."""

from __future__ import annotations

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
    assert result.returncode == 0, (
        f"snakemake dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )


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
    assert result.returncode == 0, (
        f"snakemake mock run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )
    assert (work / "reports" / "cycle_01.md").exists()
