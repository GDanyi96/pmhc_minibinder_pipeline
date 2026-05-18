"""Stage 2 designs orchestrator -- splice + MPNN + AF2 fan-in + halt gate.

Mirrors `scripts/run_stage1.py`:
- Reads Stage 1's stage1_summary.json (verdict must be PASS).
- Splices each Stage 1 backbone onto the cleaned pMHC.
- Drives ProteinMPNN per spliced PDB (or copies mock FASTAs).
- Aggregates and ranks by NLL, takes the top fan_in_top_n.
- Drives ColabFold on the top-N (or copies mock predictions).
- Runs compute_metrics on every prediction.
- Applies the halt gate (fraction_pass_intermediate >= 0.10).

Halt rule for cycle 2 is calibration-only: pipeline-health floor, tightened
in cycle 3 once partial diffusion lands. See specs/stage2_proteinmpnn_af2.md.
"""

from __future__ import annotations

import os

# RunPod base sets LD_LIBRARY_PATH=/usr/local/cuda/lib64 (CUDA 13); JAX and
# PyTorch both SIGSEGV cuDNN under that env. Pop before any CUDA-linked
# imports or subprocess spawns. Mirrors run_stage1.py:22 / run_controls.py:24.
os.environ.pop("LD_LIBRARY_PATH", None)

import argparse
import json
import logging
import shutil
import sys
import time
from pathlib import Path
from typing import Any, Literal

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from workflow.scripts import (
    aggregate_mpnn_outputs,
    compute_metrics,
    splice_binder,
)
from workflow.scripts.prep_target import TargetManifest

logger = logging.getLogger("run_stage2")
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

_HALT_THRESHOLD = 0.10
_TIGHT_IPAE_MAX = 7.0
_TIGHT_IPLDDT_MIN = 92.0

MOCK_FIXTURES_DIR = Path("tests/fixtures/stage2/designs")
DEFAULT_MPNN_CONFIG = Path("configs/proteinmpnn_default.yaml")
DEFAULT_AF2_CONFIG = Path("configs/af2_stage2.yaml")
DEFAULT_SEEDS_CONFIG = Path("configs/seeds.yaml")

HaltStatus = Literal["pass", "halt"]


def enforce_halt_gate(summary: dict[str, object]) -> tuple[HaltStatus, str]:
    """Single-rule gate on fraction_pass_intermediate.

    Mirrors `scripts.run_stage1.enforce_halt_gate`.
    """
    halt_rule = summary.get("halt_rule")
    if not isinstance(halt_rule, dict):
        return "halt", "HALT: stage2_summary.json missing halt_rule block"
    verdict = halt_rule.get("verdict")
    if verdict == "PASS":
        return "pass", ""
    observed = halt_rule.get("observed")
    threshold = halt_rule.get("threshold")
    return "halt", (
        f"HALT: fraction_pass_intermediate {observed} < {threshold} "
        "(calibration-only threshold for cycle 2; see specs/stage2_proteinmpnn_af2.md)"
    )


def _load_target(manifest_path: Path) -> TargetManifest:
    raw = yaml.safe_load(manifest_path.read_text())
    return TargetManifest.model_validate(raw)


def _resolve_cleaned_pdb(cleaned_pdb: str, manifest_path: Path) -> Path:
    candidate = Path(cleaned_pdb)
    if candidate.is_absolute():
        return candidate
    next_to_manifest = manifest_path.parent / candidate
    if next_to_manifest.exists():
        return next_to_manifest.resolve()
    return candidate


def _splice_all(backbone_pdbs: list[Path], cleaned_pmhc: Path, spliced_dir: Path) -> list[Path]:
    spliced_dir.mkdir(parents=True, exist_ok=True)
    spliced_paths: list[Path] = []
    for backbone in backbone_pdbs:
        out_pdb = spliced_dir / backbone.name
        splice_binder.splice_binder_onto_pmhc(backbone, cleaned_pmhc, out_pdb)
        spliced_paths.append(out_pdb)
    return spliced_paths


def _copy_mock_mpnn_outputs(backbones: list[Path], per_design_root: Path) -> list[Path]:
    src_dir = MOCK_FIXTURES_DIR / "mpnn_outputs"
    fastas: list[Path] = []
    for backbone in backbones:
        dest_dir = per_design_root / backbone.stem
        dest_dir.mkdir(parents=True, exist_ok=True)
        src = src_dir / f"{backbone.stem}.fa"
        dest = dest_dir / f"{backbone.stem}.fa"
        shutil.copy(src, dest)
        fastas.append(dest)
    return fastas


def _copy_mock_af2_predictions(top_records: list[dict[str, Any]], predictions_root: Path) -> None:
    src_root = MOCK_FIXTURES_DIR / "af2_predictions"
    predictions_root.mkdir(parents=True, exist_ok=True)
    for rec in top_records:
        pred_id = f"{rec['design_id']}_seq{rec['seq_id']}"
        src = src_root / pred_id
        dest = predictions_root / pred_id
        if dest.exists():
            shutil.rmtree(dest)
        shutil.copytree(src, dest)


def _compute_metrics_for_predictions(
    top_records: list[dict[str, Any]], predictions_root: Path
) -> list[dict[str, Any]]:
    metric_records: list[dict[str, Any]] = []
    for rec in top_records:
        pred_id = f"{rec['design_id']}_seq{rec['seq_id']}"
        pred_dir = predictions_root / pred_id
        m = compute_metrics.metrics_for_control(pred_dir, None)
        metric_records.append(
            {
                "design_id": rec["design_id"],
                "seq_id": rec["seq_id"],
                "iPAE": m["ipae"],
                "ipLDDT": m["iplddt"],
                "BSA": m["bsa"],
            }
        )
    return metric_records


def _write_jsonl(records: list[dict[str, Any]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        for r in records:
            fh.write(json.dumps(r, sort_keys=True) + "\n")


def _build_summary(
    cycle: int,
    target_id: str,
    n_backbones: int,
    n_sequences: int,
    metric_records: list[dict[str, Any]],
    af2_cfg: dict[str, Any],
    started_at: float,
    mock: bool,
) -> dict[str, Any]:
    cut_ipae = af2_cfg["halt_cut_ipae_max"]
    cut_iplddt = af2_cfg["halt_cut_iplddt_min"]
    threshold = af2_cfg["halt_threshold_fraction"]
    n_total = len(metric_records)
    n_pass = sum(1 for r in metric_records if r["iPAE"] < cut_ipae and r["ipLDDT"] > cut_iplddt)
    n_pass_tight = sum(
        1 for r in metric_records if r["iPAE"] < _TIGHT_IPAE_MAX and r["ipLDDT"] > _TIGHT_IPLDDT_MIN
    )
    observed = n_pass / n_total if n_total else 0.0
    fraction_tight = n_pass_tight / n_total if n_total else 0.0
    verdict = "PASS" if observed >= threshold else "FAIL"
    top_10 = sorted(metric_records, key=lambda r: r["iPAE"])[:10]
    return {
        "cycle": cycle,
        "target_id": target_id,
        "n_backbones_in": n_backbones,
        "n_sequences_designed": n_sequences,
        "n_af2_folded": n_total,
        "n_pass_intermediate": n_pass,
        "fraction_pass_intermediate": observed,
        "fraction_pass_tight": fraction_tight,
        "halt_rule": {
            "name": "fraction_pass_intermediate",
            "threshold": threshold,
            "observed": observed,
            "verdict": verdict,
            "note": "Threshold is calibration-only for cycle 2; tighten in cycle 3",
        },
        "top_10_by_ipae": top_10,
        "wall_minutes_total": round((time.time() - started_at) / 60.0, 4),
        "mock": mock,
    }


def run(
    stage1_summary_path: Path,
    target_manifest: Path,
    mpnn_config: Path,
    af2_config: Path,
    seeds_config: Path,
    cycle: int,
    out_dir: Path,
    mock: bool,
    metrics_only: bool,
) -> int:
    started_at = time.time()

    if metrics_only:
        summary_path = out_dir / "stage2_summary.json"
        if not summary_path.exists():
            logger.error("--metrics-only: %s not found", summary_path)
            return 1
        summary = json.loads(summary_path.read_text())
        status, message = enforce_halt_gate(summary)
        if status == "halt":
            logger.error(message)
            return 1
        return 0

    if not stage1_summary_path.exists():
        logger.error("stage1 summary missing: %s", stage1_summary_path)
        return 1
    stage1 = json.loads(stage1_summary_path.read_text())
    if stage1.get("halt_rule", {}).get("verdict") != "PASS":
        logger.error("Stage 1 verdict is not PASS at %s; cannot proceed", stage1_summary_path)
        return 1

    target = _load_target(target_manifest)
    cleaned_pmhc = _resolve_cleaned_pdb(target.primary.cleaned_pdb, target_manifest)
    af2_cfg = yaml.safe_load(af2_config.read_text())
    seeds = aggregate_mpnn_outputs._load_seeds(seeds_config)

    # Locate Stage 1 backbones. In mock mode the fixture ships its own
    # stage1_backbones/ directory next to stage1_summary.json; in real mode
    # backbones live under <stage1_summary_dir>/designs/.
    if mock:
        backbone_dir = stage1_summary_path.parent / "stage1_backbones"
    else:
        backbone_dir = stage1_summary_path.parent / "designs"
    backbones = sorted(backbone_dir.glob("design_*.pdb"))
    if not backbones:
        logger.error("no stage 1 backbones found under %s", backbone_dir)
        return 1

    out_dir.mkdir(parents=True, exist_ok=True)
    spliced_dir = out_dir / "proteinmpnn_designs" / "spliced"
    spliced_pdbs = _splice_all(backbones, cleaned_pmhc, spliced_dir)
    (spliced_dir / ".done").touch()

    per_design_root = out_dir / "proteinmpnn_designs" / "per_design"
    if mock:
        per_design_fastas = _copy_mock_mpnn_outputs(backbones, per_design_root)
    else:
        msg = "real ProteinMPNN invocation is the pod's job; this PR ships mock-CI-green only"
        raise NotImplementedError(msg)

    design_records: list[dict[str, str | int]] = [
        {
            "design_id": backbone.stem,
            "backbone_pdb": str(backbone),
            "spliced_pdb": str(spliced),
        }
        for backbone, spliced in zip(backbones, spliced_pdbs, strict=True)
    ]
    sequences_jsonl = out_dir / "proteinmpnn_designs" / "sequences.jsonl"
    aggregate_mpnn_outputs.aggregate(
        per_design_fastas=per_design_fastas,
        design_records=design_records,
        cycle=cycle,
        seeds=seeds,
        out_jsonl=sequences_jsonl,
    )

    records = [json.loads(line) for line in sequences_jsonl.read_text().splitlines()]
    records.sort(key=lambda r: r["mpnn_nll"])
    n_total = len(records)
    fan_in_top_n = min(n_total, 40) if mock else int(af2_cfg["fan_in_top_n"])
    top_records = records[:fan_in_top_n]

    predictions_root = out_dir / "af2_designs" / "predictions"
    if mock:
        _copy_mock_af2_predictions(top_records, predictions_root)
    else:
        msg = "real ColabFold invocation is the pod's job; this PR ships mock-CI-green only"
        raise NotImplementedError(msg)

    metric_records = _compute_metrics_for_predictions(top_records, predictions_root)
    _write_jsonl(metric_records, out_dir / "af2_designs" / "metrics.jsonl")

    summary = _build_summary(
        cycle=cycle,
        target_id=target.primary.target_id,
        n_backbones=len(backbones),
        n_sequences=n_total,
        metric_records=metric_records,
        af2_cfg=af2_cfg,
        started_at=started_at,
        mock=mock,
    )
    summary_path = out_dir / "stage2_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    status, message = enforce_halt_gate(summary)
    if status == "halt":
        logger.error(message)
        return 1
    logger.info(
        "Stage 2 halt gate passed: %d/%d AF2 folds at intermediate cut (%.3f >= %.2f)",
        summary["n_pass_intermediate"],
        summary["n_af2_folded"],
        summary["fraction_pass_intermediate"],
        _HALT_THRESHOLD,
    )
    return 0


def _default_stage1_summary(cycle: int, mock: bool) -> Path:
    if mock:
        return MOCK_FIXTURES_DIR / "stage1_summary.json"
    return Path(f"results/cycle_{cycle:02d}/stage1/rfdiffusion/stage1_summary.json")


def _default_target_manifest(cycle: int, mock: bool) -> Path:
    if mock:
        return MOCK_FIXTURES_DIR / "mock_target.yaml"
    return Path(f"results/cycle_{cycle:02d}/stage0/target.yaml")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cycle", type=int, default=2)
    parser.add_argument("--stage1-summary", type=Path, default=None)
    parser.add_argument("--target-manifest", type=Path, default=None)
    parser.add_argument("--mpnn-config", type=Path, default=DEFAULT_MPNN_CONFIG)
    parser.add_argument("--af2-config", type=Path, default=DEFAULT_AF2_CONFIG)
    parser.add_argument("--seeds-config", type=Path, default=DEFAULT_SEEDS_CONFIG)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--mock", action="store_true")
    parser.add_argument(
        "--metrics-only",
        action="store_true",
        help="Skip splice + MPNN + ColabFold + compute_metrics; re-read existing stage2_summary.json.",
    )
    args = parser.parse_args(argv)

    stage1_summary = args.stage1_summary or _default_stage1_summary(args.cycle, args.mock)
    target_manifest = args.target_manifest or _default_target_manifest(args.cycle, args.mock)
    out_dir = args.out_dir or Path(f"results/cycle_{args.cycle:02d}/stage2")

    return run(
        stage1_summary_path=stage1_summary,
        target_manifest=target_manifest,
        mpnn_config=args.mpnn_config,
        af2_config=args.af2_config,
        seeds_config=args.seeds_config,
        cycle=args.cycle,
        out_dir=out_dir,
        mock=args.mock,
        metrics_only=args.metrics_only,
    )


if __name__ == "__main__":
    sys.exit(main())
