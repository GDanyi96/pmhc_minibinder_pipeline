# mypy: disable-error-code="no-untyped-call,attr-defined,no-any-return,misc"
"""Cycle 03 Stage 1 merge -- combine sub-runs A + B into the canonical output.

Sub-runs A (BAKER scaffolds) and B (design_2079 seed) each emit their own
``designs/`` + ``designs.jsonl``. This step concatenates them into the single
handoff that Stage 2 already consumes unchanged:
``results/cycle_NN/stage1/rfdiffusion/{designs/, designs.jsonl, stage1_summary.json}``.
The geometry-pass halt gate runs once, here, on the merged set (mirrors
scripts/run_stage1.py semantics: PASS when fraction_geometry_pass >= 0.10,
calibration-only). Placement quality is enforced separately downstream by
the cycle-03 contact filter (workflow/scripts/contact_filter.py).
"""

from __future__ import annotations

import os

os.environ.pop("LD_LIBRARY_PATH", None)

import argparse
import json
import logging
import shutil
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.run_stage1 import HALT_THRESHOLD, enforce_halt_gate
from workflow.scripts.run_rfdiffusion import _length_dist, _load_seeds

logger = logging.getLogger("merge_stage1_subruns")
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

DEFAULT_SEEDS_YAML = Path("configs/seeds.yaml")


def _read_designs(subrun_dir: Path) -> list[dict[str, Any]]:
    jsonl = subrun_dir / "designs.jsonl"
    if not jsonl.exists():
        raise FileNotFoundError(f"sub-run designs manifest missing: {jsonl}")
    return [json.loads(line) for line in jsonl.read_text().splitlines() if line.strip()]


def run(
    subrun_dirs: list[Path],
    seeds_yaml: Path,
    cycle: int,
    out_dir: Path,
    mock: bool,
) -> int:
    started_at = time.time()
    designs_dir = out_dir / "designs"
    designs_dir.mkdir(parents=True, exist_ok=True)

    merged: list[dict[str, Any]] = []
    target_id = "unknown"
    for subrun_dir in subrun_dirs:
        for rec in _read_designs(subrun_dir):
            src = Path(rec["pdb_path"])
            if src.exists():
                dst = designs_dir / src.name
                if not dst.exists():
                    shutil.copy(src, dst)
                rec = {**rec, "pdb_path": str(dst)}
            target_id = rec.get("target_id", target_id)
            merged.append(rec)

    out_jsonl = out_dir / "designs.jsonl"
    with out_jsonl.open("w") as fh:
        for rec in merged:
            fh.write(json.dumps(rec, sort_keys=True) + "\n")

    n_completed = len(merged)
    n_pass = sum(1 for r in merged if r.get("geometry_pass"))
    fraction = (n_pass / n_completed) if n_completed else 0.0
    verdict = "PASS" if fraction >= HALT_THRESHOLD else "FAIL"

    seeds = _load_seeds(seeds_yaml)
    seed_lo, seed_hi = seeds.reserved[f"cycle_{cycle:02d}_rfdiffusion"]
    summary: dict[str, Any] = {
        "cycle": cycle,
        "target_id": target_id,
        "subruns": [d.name for d in subrun_dirs],
        "n_attempted": n_completed,
        "n_completed": n_completed,
        "n_geometry_pass": n_pass,
        "fraction_geometry_pass": fraction,
        "length_dist": _length_dist([int(r["binder_length"]) for r in merged]),
        "seed_range": [seed_lo, seed_hi],
        "halt_rule": {
            "name": "geometry_pass_rate",
            "threshold": HALT_THRESHOLD,
            "observed": fraction,
            "verdict": verdict,
            "note": (
                "Cycle 03 merged geometry gate (calibration-only); placement "
                "quality enforced downstream by contact_filter.py"
            ),
        },
        "wall_minutes_total": round((time.time() - started_at) / 60.0, 4),
        "mock": mock,
    }
    summary_path = out_dir / "stage1_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    logger.info(
        "merged %d designs from %d sub-runs -> %s (verdict=%s)",
        n_completed,
        len(subrun_dirs),
        summary_path,
        verdict,
    )

    status, message = enforce_halt_gate(summary)
    if status == "halt":
        logger.error(message)
        return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cycle", type=int, default=3)
    parser.add_argument("--subrun-dir", type=Path, action="append", default=None)
    parser.add_argument("--seeds-yaml", type=Path, default=DEFAULT_SEEDS_YAML)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--mock", action="store_true")
    args = parser.parse_args(argv)

    stage1_dir = Path(f"results/cycle_{args.cycle:02d}/stage1")
    subrun_dirs = args.subrun_dir or [stage1_dir / "subrun_a", stage1_dir / "subrun_b"]
    out_dir = args.out_dir or (stage1_dir / "rfdiffusion")

    return run(
        subrun_dirs=subrun_dirs,
        seeds_yaml=args.seeds_yaml,
        cycle=args.cycle,
        out_dir=out_dir,
        mock=args.mock,
    )


if __name__ == "__main__":
    sys.exit(main())
