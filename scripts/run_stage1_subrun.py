# mypy: disable-error-code="no-untyped-call,attr-defined,no-any-return,misc"
"""Cycle 03 Stage 1 sub-run orchestrator (partial diffusion).

One script drives both partial-diffusion sub-runs (parameterized by
``--subrun``):

* ``a`` -- BAKER scaffold library. Align each scaffold onto our reference
  (align_scaffolds.py), then partial-diffuse one binder backbone per
  scaffold.
* ``b`` -- our cycle-02 hero seed complex. Partial-diffuse ``num_designs``
  binder backbones from the single seed.

Each sub-run emits, under ``out_dir`` (``results/cycle_NN/stage1/subrun_X/``):
- ``designs/design_{seed}.pdb`` -- 4-chain complexes (motif A/B/C + binder D)
- ``designs.jsonl`` -- one geometry-scored record per design
- ``subrun_summary.json`` -- counts + geometry-pass fraction (no halt
  verdict; the gate runs once on the merged set, merge_stage1_subruns.py).

Seeds come from the cycle's reserved rfdiffusion range; the per-subrun
``seed_offset`` (configs/rfdiffusion_subrun_{a,b}.yaml) keeps A and B in
disjoint slots. Reuses the Stage 1 geometry/seed helpers in
workflow.scripts.run_rfdiffusion so the design contract is identical.
"""

from __future__ import annotations

import os

os.environ.pop("LD_LIBRARY_PATH", None)

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from workflow.scripts import align_scaffolds, partial_diffuse
from workflow.scripts.run_rfdiffusion import (
    RFDIFFUSION_DIR,
    _compute_seed,
    _emit_designs_jsonl,
    _file_sha256,
    _geometry_pass,
    _git_sha,
    _hotspot_ca_xyz,
    _hotspot_tokens,
    _length_dist,
    _load_seeds,
    _load_target,
)

logger = logging.getLogger("run_stage1_subrun")
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

DEFAULT_SEEDS_YAML = Path("configs/seeds.yaml")
MOCK_TARGET_MANIFEST = Path("tests/fixtures/stage1/mock_target.yaml")
MOCK_REFERENCE_PDB = Path("tests/fixtures/stage1/mock_clean.pdb")
MOCK_SCAFFOLD_GLOB = "tests/fixtures/baker_library_mock/scaf*.pdb"
MOCK_SEED_PDB = Path("tests/fixtures/design_2079_mock_seed.pdb")
MOCK_BINDER_LENGTH = 80
_MOCK_SUBRUN_B_DESIGNS = 2


def _resolve_inputs(
    subrun: str,
    cfg: dict[str, Any],
    out_dir: Path,
    reference_pdb: Path,
    mock: bool,
) -> list[Path]:
    """Return the per-design input PDB list (aligned scaffolds for A, the
    repeated seed complex for B)."""
    if subrun == "a":
        scaffold_glob = MOCK_SCAFFOLD_GLOB if mock else (cfg.get("scaffolds") or {})["aligned_glob"]
        # Real sub-run A uses BAKER's fused-chain library, aligned against the
        # truncated reference (chain B = HLA, C = peptide) and rewritten to our
        # A=binder / B=HLA / C=peptide layout. Mock keeps the lightweight
        # toy-scaffold aligner + mock_clean.pdb so the cycle-99 DAG stays green.
        aligned = align_scaffolds.align_scaffolds(
            str(scaffold_glob), reference_pdb, out_dir / "aligned", baker_layout=not mock
        )
        n = len(aligned) if mock else min(int(cfg["inference"]["num_designs"]), len(aligned))
        return aligned[:n]
    seed_pdb = MOCK_SEED_PDB if mock else Path(cfg["inference"]["input_pdb"])
    n = _MOCK_SUBRUN_B_DESIGNS if mock else int(cfg["inference"]["num_designs"])
    return [seed_pdb] * n


def run(
    subrun: str,
    config_yaml: Path,
    seeds_yaml: Path,
    target_manifest: Path,
    reference_pdb: Path,
    cycle: int,
    out_dir: Path,
    mock: bool,
) -> int:
    started_at = time.time()
    cfg = yaml.safe_load(config_yaml.read_text())
    seeds = _load_seeds(seeds_yaml)
    seed_offset = int(cfg.get("seed_offset", 0))
    partial_T = int(cfg.get("partial_T", partial_diffuse.DEFAULT_PARTIAL_T))
    noise_scale_ca = int(cfg.get("denoiser", {}).get("noise_scale_ca", 0))
    ckpt_path = Path(cfg["inference"]["ckpt_override_path"])

    manifest = _load_target(target_manifest)
    target = manifest.primary
    cleaned_pdb = Path(target.cleaned_pdb)
    hotspot_tokens = _hotspot_tokens(target)
    hotspot_xyz = _hotspot_ca_xyz(target, cleaned_pdb)
    fixed_chains = frozenset(cid for cid, chain in target.chains.items() if chain.role != "binder")
    length_range = (target.binder.length_min, target.binder.length_max)

    designs_dir = out_dir / "designs"
    designs_dir.mkdir(parents=True, exist_ok=True)
    hydra_dir = out_dir / "hydra_logs"

    inputs = _resolve_inputs(subrun, cfg, out_dir, reference_pdb, mock)
    if not inputs:
        logger.error("sub-run %s: no input PDBs resolved", subrun)
        return 1

    design_paths: list[Path] = []
    for i, input_pdb in enumerate(inputs):
        seed = _compute_seed(cycle, seed_offset + i, seeds)
        out_pdb = designs_dir / f"design_{seed}.pdb"
        if not out_pdb.exists():
            if mock:
                partial_diffuse._synthesize_mock_design(out_pdb, reference_pdb, MOCK_BINDER_LENGTH)
            else:
                partial_diffuse.partial_diffuse_one(
                    input_pdb=input_pdb,
                    out_prefix=out_pdb.with_suffix(""),
                    seed=seed,
                    ckpt=ckpt_path,
                    partial_T=partial_T,
                    noise_scale_ca=noise_scale_ca,
                    hydra_dir=hydra_dir / f"design_{seed}",
                )
        design_paths.append(out_pdb)

    ckpt_sha = _file_sha256(ckpt_path)
    rfdiff_sha = _git_sha(RFDIFFUSION_DIR)

    records: list[dict[str, Any]] = []
    for i, pdb in enumerate(design_paths):
        if not pdb.exists():
            continue
        geom = _geometry_pass(pdb, hotspot_xyz, length_range, fixed_chains)
        seed = _compute_seed(cycle, seed_offset + i, seeds)
        records.append(
            {
                "design_id": pdb.stem,
                "pdb_path": str(pdb),
                "target_id": target.target_id,
                "subrun": subrun,
                "binder_length": geom["binder_length"],
                "hotspot_residues": hotspot_tokens,
                "contigmap": None,
                "partial_T": partial_T,
                "noise_scale_ca": float(noise_scale_ca),
                "noise_scale_frame": 0.0,
                "checkpoint": ckpt_path.name,
                "checkpoint_sha256": ckpt_sha,
                "rfdiffusion_git_sha": rfdiff_sha,
                "seed": seed,
                "ca_contact_to_hotspots_n": geom["ca_contact_to_hotspots_n"],
                "geometry_pass": geom["geometry_pass"],
                "fail_reasons": geom["fail_reasons"],
            }
        )

    _emit_designs_jsonl(records, out_dir / "designs.jsonl")

    n_pass = sum(1 for r in records if r["geometry_pass"])
    n_completed = len(records)
    fraction = (n_pass / n_completed) if n_completed else 0.0
    summary = {
        "cycle": cycle,
        "subrun": subrun,
        "target_id": target.target_id,
        "n_attempted": len(inputs),
        "n_completed": n_completed,
        "n_geometry_pass": n_pass,
        "fraction_geometry_pass": fraction,
        "length_dist": _length_dist([r["binder_length"] for r in records]),
        "partial_T": partial_T,
        "seed_offset": seed_offset,
        "rfdiffusion_git_sha": rfdiff_sha,
        "weights_sha256": {ckpt_path.name: ckpt_sha},
        "wall_minutes_total": round((time.time() - started_at) / 60.0, 4),
        "mock": mock,
    }
    summary_path = out_dir / "subrun_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    logger.info(
        "sub-run %s: %d/%d geometry-pass -> %s",
        subrun,
        n_pass,
        n_completed,
        summary_path,
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subrun", choices=["a", "b"], required=True)
    parser.add_argument("--cycle", type=int, default=3)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--seeds-yaml", type=Path, default=DEFAULT_SEEDS_YAML)
    parser.add_argument("--target-manifest", type=Path, default=None)
    parser.add_argument("--reference-pdb", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--mock", action="store_true")
    args = parser.parse_args(argv)

    config_yaml = args.config or Path(f"configs/rfdiffusion_subrun_{args.subrun}.yaml")
    target_manifest = args.target_manifest or (
        MOCK_TARGET_MANIFEST
        if args.mock
        else Path(f"results/cycle_{args.cycle:02d}/stage0/target.yaml")
    )
    reference_pdb = args.reference_pdb or (
        MOCK_REFERENCE_PDB if args.mock else Path("data/targets/3hpj_clean.pdb")
    )
    out_dir = args.out_dir or Path(f"results/cycle_{args.cycle:02d}/stage1/subrun_{args.subrun}")

    return run(
        subrun=args.subrun,
        config_yaml=config_yaml,
        seeds_yaml=args.seeds_yaml,
        target_manifest=target_manifest,
        reference_pdb=reference_pdb,
        cycle=args.cycle,
        out_dir=out_dir,
        mock=args.mock,
    )


if __name__ == "__main__":
    sys.exit(main())
