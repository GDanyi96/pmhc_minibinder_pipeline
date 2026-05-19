"""Stage 1 orchestrator -- RFdiffusion native run + halt-gate enforcement.

Mirrors the shape of `scripts/run_controls.py`:
- Loads the Stage-0 cycle-scoped TargetManifest.
- Calls `workflow.scripts.run_rfdiffusion.run_rfdiffusion`.
- Reads the emitted stage1_summary.json.
- Translates the halt verdict into the process exit code (PASS=0, FAIL=1).

Halt rule for cycle 2 is calibration-only: fraction_geometry_pass >= 0.50.
Per BAKER_LAB_2025 + HADRUP_JENKINS_2025 we have zero empirical data on
the per-contigmap pass rate, so 0.50 is a "not catastrophically broken"
floor; tightened in cycle 3 once cycle-2 distribution is known.
"""

from __future__ import annotations

import os

# RunPod base sets LD_LIBRARY_PATH=/usr/local/cuda/lib64 (CUDA 13), which
# overrides PyTorch's bundled CUDA libs and SIGSEGVs cuDNN in the
# RFdiffusion subprocess we spawn.
os.environ.pop("LD_LIBRARY_PATH", None)

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Literal

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from workflow.scripts.run_rfdiffusion import run_rfdiffusion

logger = logging.getLogger("run_stage1")
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

HALT_THRESHOLD = 0.50
MOCK_FIXTURES_DIR = Path("tests/fixtures/stage1")
MOCK_TARGET_MANIFEST = MOCK_FIXTURES_DIR / "mock_target.yaml"
DEFAULT_RFDIFF_YAML = Path("configs/rfdiffusion_default.yaml")
DEFAULT_SEEDS_YAML = Path("configs/seeds.yaml")

HaltStatus = Literal["pass", "halt"]


def enforce_halt_gate(summary: dict[str, object]) -> tuple[HaltStatus, str]:
    """Single-rule gate on the geometry-pass fraction.

    Mirrors the structure of `scripts.run_controls.enforce_halt_gate` so
    the orchestrator pattern stays consistent across stages.
    """
    halt_rule = summary.get("halt_rule")
    if not isinstance(halt_rule, dict):
        return "halt", "HALT: stage1_summary.json missing halt_rule block"
    verdict = halt_rule.get("verdict")
    if verdict == "PASS":
        return "pass", ""
    observed = halt_rule.get("observed")
    threshold = halt_rule.get("threshold")
    return "halt", (
        f"HALT: fraction_geometry_pass {observed} < {threshold} "
        "(calibration-only threshold for cycle 2; see specs/stage1_rfdiffusion.md)"
    )


def _default_target_manifest(cycle: int, mock: bool) -> Path:
    if mock:
        return MOCK_TARGET_MANIFEST
    return Path(f"results/cycle_{cycle:02d}/stage0/target.yaml")


def _ensure_mock_inputs(target_manifest: Path) -> None:
    """In mock mode, make sure the fixture manifest exists. We don't copy
    it anywhere -- run_rfdiffusion reads it directly from
    tests/fixtures/stage1/ and references mock_clean.pdb relative to that
    directory."""
    if not target_manifest.exists():
        raise FileNotFoundError(
            f"mock target manifest missing at {target_manifest} -- "
            "run `uv run python tests/fixtures/stage1/_make_fixtures.py` first"
        )


def run(
    target_manifest: Path,
    rfdiff_yaml: Path,
    seeds_yaml: Path,
    cycle: int,
    out_dir: Path,
    mock: bool,
    skip_subprocess: bool,
) -> int:
    if mock:
        _ensure_mock_inputs(target_manifest)

    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = run_rfdiffusion(
        target_manifest=target_manifest,
        rfdiff_yaml=rfdiff_yaml,
        seeds_yaml=seeds_yaml,
        cycle=cycle,
        out_dir=out_dir,
        mock=mock,
        halt_threshold=HALT_THRESHOLD,
        mock_fixtures_dir=MOCK_FIXTURES_DIR if mock else None,
        skip_subprocess=skip_subprocess,
    )

    summary = json.loads(summary_path.read_text())
    status, message = enforce_halt_gate(summary)
    if status == "halt":
        logger.error(message)
        return 1
    logger.info(
        "Stage 1 halt gate passed: %s/%s designs geometry-pass (%.3f >= %.2f)",
        summary["n_geometry_pass"],
        summary["n_completed"],
        summary["fraction_geometry_pass"],
        HALT_THRESHOLD,
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cycle", type=int, default=2)
    parser.add_argument("--target-manifest", type=Path, default=None)
    parser.add_argument("--rfdiff-yaml", type=Path, default=DEFAULT_RFDIFF_YAML)
    parser.add_argument("--seeds-yaml", type=Path, default=DEFAULT_SEEDS_YAML)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--mock", action="store_true")
    parser.add_argument(
        "--skip-subprocess",
        action="store_true",
        help=(
            "Skip the RFdiffusion subprocess (and mock-fixture copy); "
            "re-enumerate existing design_*.pdb files in out-dir/designs "
            "and rewrite stage1_summary.json. Recovery path for cycle 02 "
            "writer/reader filename drift (see docs/known_traps.md)."
        ),
    )
    args = parser.parse_args(argv)

    target_manifest = args.target_manifest or _default_target_manifest(args.cycle, args.mock)
    out_dir = args.out_dir or Path(f"results/cycle_{args.cycle:02d}/stage1/rfdiffusion")

    return run(
        target_manifest=target_manifest,
        rfdiff_yaml=args.rfdiff_yaml,
        seeds_yaml=args.seeds_yaml,
        cycle=args.cycle,
        out_dir=out_dir,
        mock=args.mock,
        skip_subprocess=args.skip_subprocess,
    )


if __name__ == "__main__":
    sys.exit(main())
