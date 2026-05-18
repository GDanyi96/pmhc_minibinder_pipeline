# mypy: disable-error-code="no-untyped-call,attr-defined,no-any-return,misc"
"""Stage 1 -- Native RFdiffusion subprocess driver.

Generates de novo binder backbones against the cleaned pMHC-I target using
RFdiffusion's complex binder mode (HADRUP_JENKINS_2025, BAKER_LAB_2025).

Subprocess pattern is governed by ``_SEED_THREADING_MODE``:

    "per_design"        -- spawn one ``run_inference.py`` subprocess per
                          design (200 calls x ~10-20 s model-load
                          overhead). Safe default until the pod recon
                          below confirms threading.
    "single_subprocess" -- spawn one subprocess for the whole batch with
                          ``inference.num_designs=N inference.random_seed=base``,
                          relying on RFdiffusion to seed designs
                          ``base, base+1, ..., base+N-1``. Saves ~30-60 min
                          of model-load overhead per cycle.

Pod recon to decide the mode (zero-LOC; ~5 min on the pod):

    uv run python /workspace/RFdiffusion/scripts/run_inference.py \\
        inference.output_prefix=/tmp/seedtest/design \\
        inference.num_designs=3 inference.random_seed=2000 \\
        inference.ckpt_override_path=/workspace/models/RFdiffusion/Complex_base_ckpt.pt \\
        contigmap.contigs=[\\"A1-275/0 B1-99/0 C1-9/0 80-80\\"] \\
        ppi.hotspot_res=[A65,A150,C4]
    grep -i seed /tmp/seedtest/design_*.trb 2>&1 | head -20

If the three designs show distinct seeds 2000/2001/2002, flip
``_SEED_THREADING_MODE`` to ``"single_subprocess"`` and document in
``docs/known_traps.md``. Otherwise keep ``"per_design"`` as the proven
fallback.

This module ships only the helper machinery and the public
``run_rfdiffusion`` entry point. The orchestrator (``scripts/run_stage1.py``)
owns argparsing, halt-gate enforcement, and JSON emission glue.
"""

from __future__ import annotations

import os

# RunPod base sets LD_LIBRARY_PATH=/usr/local/cuda/lib64 (CUDA 13), which
# overrides PyTorch's bundled CUDA libs and SIGSEGVs cuDNN. Any subprocess
# we spawn (RFdiffusion's run_inference.py) inherits this process's env,
# so popping here propagates the unset to the child.
os.environ.pop("LD_LIBRARY_PATH", None)

import argparse
import json
import logging
import subprocess
import sys
import time
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Literal

import numpy as np
import yaml
from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
from pydantic import BaseModel

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from workflow.scripts.prep_target import (
    ResolvedHotspot,
    ResolvedTarget,
    TargetManifest,
)

logger = logging.getLogger("run_rfdiffusion")
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

_SEED_THREADING_MODE: Literal["single_subprocess", "per_design"] = "per_design"

CA_HOTSPOT_CONTACT_DISTANCE = 10.0
CA_INTERNAL_CLASH_DISTANCE = 3.5
LENGTH_TOLERANCE = 5
MIN_HOTSPOT_CONTACTS = 3

RFDIFFUSION_DIR = Path(os.environ.get("RFDIFFUSION_DIR", "/workspace/RFdiffusion"))
RFDIFFUSION_MODELS_DIR = Path(
    os.environ.get("RFDIFFUSION_MODELS_DIR", "/workspace/models/RFdiffusion")
)


def _resolve_rfdiffusion_python() -> str:
    """Return the SE3nv conda interpreter path for RFdiffusion inference.

    PROJECT_STATE.md §6 rule 4: pipeline scripts use uv, RFdiffusion uses
    SE3nv. RFdiffusion's deps (omegaconf, hydra, e3nn, dgl, SE3Transformer,
    torch 1.9+cu111) are not in the uv venv, so launching via `uv run`
    dies with ModuleNotFoundError.
    """
    path = os.environ.get(
        "RFDIFFUSION_PYTHON",
        "/workspace/miniconda3/envs/SE3nv/bin/python",
    )
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"RFDIFFUSION_PYTHON not found: {path}. "
            "Install SE3nv conda env and overlay torch==1.9.0+cu111 "
            "(see PROJECT_STATE.md §9 steps 3-5)."
        )
    return path


class SeedsConfig(BaseModel):
    formulas: dict[str, str]
    reserved: dict[str, list[int]]


def _resolve_cleaned_pdb(raw_path: str, manifest_path: Path) -> str:
    """Resolve TargetManifest.primary.cleaned_pdb to an existing file.

    Two conventions in play:
    * Mock fixture manifests bundle the cleaned PDB next to themselves;
      cleaned_pdb is a bare filename like "mock_clean.pdb".
    * Real manifests written by prep_target.py use cwd-relative paths like
      "data/targets/3hpj_clean.pdb" (resolved when snakemake runs from
      repo root).

    Strategy: if absolute, take as-is. Otherwise prefer the file next to
    the manifest, then fall back to cwd-relative. Absolute paths must
    never be baked into committed fixtures (trap #16).
    """
    candidate = Path(raw_path)
    if candidate.is_absolute():
        return str(candidate)
    next_to_manifest = manifest_path.parent / candidate
    if next_to_manifest.exists():
        return str(next_to_manifest.resolve())
    return str(candidate)


def _load_target(manifest_path: Path) -> TargetManifest:
    with manifest_path.open() as fh:
        raw = yaml.safe_load(fh)
    manifest = TargetManifest.model_validate(raw)
    manifest.primary.cleaned_pdb = _resolve_cleaned_pdb(manifest.primary.cleaned_pdb, manifest_path)
    if manifest.positive_control is not None:
        manifest.positive_control.cleaned_pdb = _resolve_cleaned_pdb(
            manifest.positive_control.cleaned_pdb, manifest_path
        )
    return manifest


def _load_seeds(seeds_yaml: Path) -> SeedsConfig:
    with seeds_yaml.open() as fh:
        return SeedsConfig.model_validate(yaml.safe_load(fh))


def _read_chain_lengths(cleaned_pdb: Path) -> dict[str, int]:
    """Return {chain_id: number_of_standard_residues} for a cleaned PDB."""
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(cleaned_pdb.stem, str(cleaned_pdb))
    model = structure[0]
    out: dict[str, int] = {}
    for chain in model.get_chains():
        residues = [r for r in chain.get_residues() if r.id[0].strip() == ""]
        out[chain.id] = len(residues)
    return out


def _build_contigmap(target: ResolvedTarget, cleaned_pdb: Path) -> str:
    """Build the RFdiffusion contigmap string from runtime chain lengths.

    Format: "A1-<hc_len>/0 B1-<b2m_len>/0 C1-<pep_len>/0 <min>-<max>".

    Chain lengths read at runtime from the cleaned PDB to avoid the cycle-1
    cleavage trap (hardcoded "A1-275 B1-99 C1-9" would silently mismatch on
    a different pMHC).
    """
    chain_lengths = _read_chain_lengths(cleaned_pdb)
    for cid in ("A", "B", "C"):
        if cid not in chain_lengths:
            raise ValueError(f"{cleaned_pdb}: chain {cid} absent")
    hc_len = chain_lengths["A"]
    b2m_len = chain_lengths["B"]
    pep_len = chain_lengths["C"]
    lo = target.binder.length_min
    hi = target.binder.length_max
    return f"A1-{hc_len}/0 B1-{b2m_len}/0 C1-{pep_len}/0 {lo}-{hi}"


def _compute_seed(cycle: int, design_index: int, seeds: SeedsConfig) -> int:
    """Apply the rfdiffusion seed formula, enforce the cycle's reserved range."""
    expected = "cycle * 1000 + design_index"
    actual = seeds.formulas.get("rfdiffusion")
    if actual != expected:
        raise ValueError(
            f"unsupported rfdiffusion seed formula {actual!r}; "
            "update run_rfdiffusion._compute_seed when adding new formulas"
        )
    seed = cycle * 1000 + design_index
    key = f"cycle_{cycle:02d}_rfdiffusion"
    if key not in seeds.reserved:
        raise ValueError(f"seeds.yaml: no reserved range for {key}")
    lo, hi = seeds.reserved[key]
    if not lo <= seed <= hi:
        raise ValueError(
            f"seed {seed} for cycle={cycle} design_index={design_index} "
            f"outside reserved range [{lo}, {hi}] in seeds.yaml"
        )
    return seed


def _binder_ca_coords(pdb_path: Path) -> np.ndarray:
    """Return (N, 3) array of Ca coords from chain A of an RFdiffusion output.

    RFdiffusion writes the binder as chain "A" inside its output PDB (per
    spec; renaming to chain D is Stage 2's splicing job). If the PDB has
    multiple chains we still take chain "A". Empty/parseless PDBs yield
    an empty (0, 3) array.
    """
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    models = list(structure.get_models())
    if not models:
        return np.empty((0, 3), dtype=np.float64)
    model = models[0]
    chain_ids = {c.id for c in model.get_chains()}
    if not chain_ids:
        return np.empty((0, 3), dtype=np.float64)
    chain_id = "A" if "A" in chain_ids else next(iter(chain_ids))
    chain = model[chain_id]
    coords: list[list[float]] = []
    for residue in chain.get_residues():
        if residue.id[0].strip() != "":
            continue
        if "CA" in residue:
            atom = residue["CA"]
            coords.append([float(c) for c in atom.get_coord()])
    if not coords:
        return np.empty((0, 3), dtype=np.float64)
    return np.asarray(coords, dtype=np.float64)


def _geometry_pass(
    pdb_path: Path,
    hotspot_ca_xyz: np.ndarray,
    length_range: tuple[int, int],
) -> dict[str, Any]:
    """Three-check geometry pass for a candidate binder backbone.

    1. >=3 binder Ca within ``CA_HOTSPOT_CONTACT_DISTANCE`` A of any hotspot.
    2. Binder length within ``[length_range[0]-5, length_range[1]+5]``.
    3. No internal binder Ca-Ca clash < ``CA_INTERNAL_CLASH_DISTANCE`` A
       (excluding sequential i, i+1 pairs -- those are the trivial backbone
       bond ~3.8 A nominal but can be slightly shorter in some outputs).
    """
    binder_xyz = _binder_ca_coords(pdb_path)
    binder_length = int(binder_xyz.shape[0])

    if binder_length == 0:
        return {
            "geometry_pass": False,
            "ca_contact_to_hotspots_n": 0,
            "binder_length": 0,
            "internal_clash_count": 0,
            "fail_reasons": ["empty_binder"],
        }

    if hotspot_ca_xyz.size > 0:
        diffs = binder_xyz[:, None, :] - hotspot_ca_xyz[None, :, :]
        per_pair_dist = np.linalg.norm(diffs, axis=-1)
        per_binder_min = per_pair_dist.min(axis=1)
        n_contacts = int((per_binder_min < CA_HOTSPOT_CONTACT_DISTANCE).sum())
    else:
        n_contacts = 0

    lo_ok = length_range[0] - LENGTH_TOLERANCE
    hi_ok = length_range[1] + LENGTH_TOLERANCE
    length_ok = lo_ok <= binder_length <= hi_ok

    clash_count = 0
    if binder_length >= 3:
        intra = np.linalg.norm(binder_xyz[:, None, :] - binder_xyz[None, :, :], axis=-1)
        n = binder_xyz.shape[0]
        mask = np.zeros((n, n), dtype=bool)
        i_idx, j_idx = np.triu_indices(n, k=2)
        mask[i_idx, j_idx] = True
        clash_count = int(((intra < CA_INTERNAL_CLASH_DISTANCE) & mask).sum())

    fail_reasons: list[str] = []
    if n_contacts < MIN_HOTSPOT_CONTACTS:
        fail_reasons.append("insufficient_hotspot_contacts")
    if not length_ok:
        fail_reasons.append("length_out_of_range")
    if clash_count > 0:
        fail_reasons.append("internal_ca_clash")

    return {
        "geometry_pass": not fail_reasons,
        "ca_contact_to_hotspots_n": n_contacts,
        "binder_length": binder_length,
        "internal_clash_count": clash_count,
        "fail_reasons": fail_reasons,
    }


def _hotspot_tokens(target: ResolvedTarget) -> list[str]:
    out: list[str] = []
    for group in ("peptide", "mhc"):
        for h in target.hotspots[group]:
            out.append(h.token)
    return out


def _hotspot_ca_xyz(target: ResolvedTarget, cleaned_pdb: Path) -> np.ndarray:
    """Read Ca coordinates of the resolved hotspot residues from the cleaned PDB."""
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(cleaned_pdb.stem, str(cleaned_pdb))
    model = structure[0]
    coords: list[list[float]] = []
    flat: list[ResolvedHotspot] = []
    for group in ("peptide", "mhc"):
        flat.extend(target.hotspots[group])
    for h in flat:
        if h.chain not in model:
            raise ValueError(f"hotspot {h.token}: chain {h.chain} absent from {cleaned_pdb}")
        chain = model[h.chain]
        key = (" ", h.resnum_clean, " ")
        if key not in chain:
            raise ValueError(
                f"hotspot {h.token}: residue {h.resnum_clean} absent from chain {h.chain}"
            )
        residue = chain[key]
        if "CA" not in residue:
            raise ValueError(f"hotspot {h.token}: residue has no Ca")
        coords.append([float(c) for c in residue["CA"].get_coord()])
    return np.asarray(coords, dtype=np.float64)


def _emit_designs_jsonl(records: Sequence[dict[str, Any]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        for record in records:
            fh.write(json.dumps(record, sort_keys=True) + "\n")


def _length_dist(lengths: Sequence[int]) -> dict[str, int | float]:
    if not lengths:
        return {"min": 0, "max": 0, "median": 0, "p25": 0, "p75": 0}
    arr = np.asarray(lengths)
    return {
        "min": int(arr.min()),
        "max": int(arr.max()),
        "median": float(np.median(arr)),
        "p25": float(np.percentile(arr, 25)),
        "p75": float(np.percentile(arr, 75)),
    }


def _git_sha(repo_dir: Path) -> str:
    if not (repo_dir / ".git").exists():
        return "unknown"
    try:
        out = subprocess.run(
            ["git", "-C", str(repo_dir), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
        return out.stdout.strip()
    except subprocess.CalledProcessError:
        return "unknown"


def _file_sha256(path: Path) -> str:
    if not path.exists():
        return "absent"
    try:
        out = subprocess.run(["sha256sum", str(path)], check=True, capture_output=True, text=True)
        return out.stdout.split()[0]
    except subprocess.CalledProcessError:
        return "unknown"


def _run_inference_per_design(
    out_dir: Path, contig: str, hotspot_tokens: list[str], seed: int, ckpt: Path
) -> None:
    rfdiffusion_python = _resolve_rfdiffusion_python()
    cmd = [
        rfdiffusion_python,
        str(RFDIFFUSION_DIR / "scripts" / "run_inference.py"),
        f"inference.output_prefix={out_dir}/design",
        "inference.num_designs=1",
        f"inference.ckpt_override_path={ckpt}",
        f"contigmap.contigs=[{contig}]",
        f"ppi.hotspot_res=[{','.join(hotspot_tokens)}]",
        "denoiser.noise_scale_ca=0",
        "denoiser.noise_scale_frame=0",
        f"inference.random_seed={seed}",
    ]
    logger.info(f"Launching RFdiffusion with interpreter: {rfdiffusion_python}")
    logger.info(f"Full cmd: {cmd}")
    subprocess.run(cmd, check=True, env={**os.environ})


def _run_inference_single_batch(
    out_dir: Path,
    contig: str,
    hotspot_tokens: list[str],
    base_seed: int,
    num_designs: int,
    ckpt: Path,
) -> None:
    rfdiffusion_python = _resolve_rfdiffusion_python()
    cmd = [
        rfdiffusion_python,
        str(RFDIFFUSION_DIR / "scripts" / "run_inference.py"),
        f"inference.output_prefix={out_dir}/design",
        f"inference.num_designs={num_designs}",
        f"inference.ckpt_override_path={ckpt}",
        f"contigmap.contigs=[{contig}]",
        f"ppi.hotspot_res=[{','.join(hotspot_tokens)}]",
        "denoiser.noise_scale_ca=0",
        "denoiser.noise_scale_frame=0",
        f"inference.random_seed={base_seed}",
    ]
    logger.info(f"Launching RFdiffusion with interpreter: {rfdiffusion_python}")
    logger.info(f"Full cmd: {cmd}")
    subprocess.run(cmd, check=True, env={**os.environ})


def _copy_mock_fixtures(out_designs_dir: Path, fixtures_dir: Path) -> int:
    """Copy mock backbone PDBs from fixtures to the per-cycle designs dir.

    Returns the number of fixtures copied. Used by mock mode in lieu of
    a real RFdiffusion subprocess.
    """
    import shutil

    out_designs_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for src in sorted(fixtures_dir.glob("mock_design_*.pdb")):
        dst_name = src.name.removeprefix("mock_")
        dst = out_designs_dir / dst_name
        if not dst.exists():
            shutil.copy(src, dst)
        n += 1
    return n


def run_rfdiffusion(
    target_manifest: Path,
    rfdiff_yaml: Path,
    seeds_yaml: Path,
    cycle: int,
    out_dir: Path,
    mock: bool,
    halt_threshold: float = 0.50,
    mock_fixtures_dir: Path | None = None,
) -> Path:
    """Run Stage 1 RFdiffusion and emit stage1_summary.json.

    Returns the path to the emitted summary. The orchestrator
    (``scripts/run_stage1.py``) is responsible for translating that into
    a process exit code via halt-gate enforcement.
    """
    started_at = time.time()

    manifest = _load_target(target_manifest)
    target = manifest.primary
    seeds = _load_seeds(seeds_yaml)

    with rfdiff_yaml.open() as fh:
        rfdiff_cfg = yaml.safe_load(fh)
    num_designs = int(rfdiff_cfg["inference"]["num_designs"])
    ckpt_path = Path(rfdiff_cfg["inference"]["ckpt_override_path"])
    partial_T = rfdiff_cfg.get("partial_T")

    cleaned_pdb = Path(target.cleaned_pdb)
    contig = _build_contigmap(target, cleaned_pdb)
    hotspot_tokens = _hotspot_tokens(target)
    hotspot_xyz = _hotspot_ca_xyz(target, cleaned_pdb)

    designs_dir = out_dir / "designs"
    designs_dir.mkdir(parents=True, exist_ok=True)

    if mock:
        fixtures_dir = mock_fixtures_dir or Path("tests/fixtures/stage1")
        n_copied = _copy_mock_fixtures(designs_dir, fixtures_dir)
        logger.info("mock: copied %d fixture PDBs from %s", n_copied, fixtures_dir)
        effective_designs = n_copied
    else:
        if _SEED_THREADING_MODE == "single_subprocess":
            base_seed = _compute_seed(cycle, 0, seeds)
            _ = _compute_seed(cycle, num_designs - 1, seeds)
            _run_inference_single_batch(
                designs_dir, contig, hotspot_tokens, base_seed, num_designs, ckpt_path
            )
        else:
            for design_index in range(num_designs):
                seed = _compute_seed(cycle, design_index, seeds)
                expected_pdb = designs_dir / f"design_{design_index:05d}.pdb"
                if expected_pdb.exists():
                    continue
                _run_inference_per_design(designs_dir, contig, hotspot_tokens, seed, ckpt_path)
        effective_designs = num_designs

    rfdiff_sha = _git_sha(RFDIFFUSION_DIR)
    ckpt_sha = _file_sha256(ckpt_path)

    records: list[dict[str, Any]] = []
    for design_index in range(effective_designs):
        pdb = designs_dir / f"design_{design_index:05d}.pdb"
        if not pdb.exists():
            continue
        geom = _geometry_pass(
            pdb, hotspot_xyz, (target.binder.length_min, target.binder.length_max)
        )
        seed = _compute_seed(cycle, design_index, seeds) if not mock else 0
        records.append(
            {
                "design_id": pdb.stem,
                "pdb_path": str(pdb),
                "target_id": target.target_id,
                "binder_length": geom["binder_length"],
                "hotspot_residues": hotspot_tokens,
                "contigmap": f"[{contig}]",
                "noise_scale_ca": 0.0,
                "noise_scale_frame": 0.0,
                "partial_T": partial_T,
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
    verdict = "PASS" if fraction >= halt_threshold else "FAIL"

    seed_lo, seed_hi = seeds.reserved[f"cycle_{cycle:02d}_rfdiffusion"]
    summary: dict[str, Any] = {
        "cycle": cycle,
        "target_id": target.target_id,
        "n_attempted": effective_designs,
        "n_completed": n_completed,
        "n_geometry_pass": n_pass,
        "fraction_geometry_pass": fraction,
        "length_dist": _length_dist([r["binder_length"] for r in records]),
        "seed_range": [seed_lo, seed_hi],
        "halt_rule": {
            "name": "geometry_pass_rate",
            "threshold": halt_threshold,
            "observed": fraction,
            "verdict": verdict,
            "note": "Threshold is calibration-only for cycle 2; tighten in cycle 3",
        },
        "rfdiffusion_git_sha": rfdiff_sha,
        "weights_sha256": {ckpt_path.name: ckpt_sha},
        "wall_minutes_total": round((time.time() - started_at) / 60.0, 2),
        "mock": mock,
    }

    summary_path = out_dir / "stage1_summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    logger.info("wrote %s -- verdict=%s observed=%.3f", summary_path, verdict, fraction)
    return summary_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-manifest", type=Path, required=True)
    parser.add_argument(
        "--rfdiff-yaml", type=Path, default=Path("configs/rfdiffusion_default.yaml")
    )
    parser.add_argument("--seeds-yaml", type=Path, default=Path("configs/seeds.yaml"))
    parser.add_argument("--cycle", type=int, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--mock", action="store_true")
    parser.add_argument("--mock-fixtures-dir", type=Path, default=None)
    args = parser.parse_args(argv)

    run_rfdiffusion(
        args.target_manifest,
        args.rfdiff_yaml,
        args.seeds_yaml,
        args.cycle,
        args.out_dir,
        args.mock,
        mock_fixtures_dir=args.mock_fixtures_dir,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
