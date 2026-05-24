# mypy: disable-error-code="no-untyped-call,attr-defined"
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
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Literal

import yaml
from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
from Bio.SeqUtils import seq1

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.run_stage1 import HALT_THRESHOLD as STAGE1_HALT_THRESHOLD
from workflow.scripts import (
    aggregate_mpnn_outputs,
    compute_metrics,
    splice_binder,
)
from workflow.scripts.prep_target import TargetManifest
from workflow.scripts.run_colabfold import _reshape_flat_outputs

logger = logging.getLogger("run_stage2")
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

_HALT_THRESHOLD = 0.10
_TIGHT_IPAE_MAX = 7.0
_TIGHT_IPLDDT_MIN = 92.0

MOCK_FIXTURES_DIR = Path("tests/fixtures/stage2/designs")
DEFAULT_MPNN_CONFIG = Path("configs/proteinmpnn_default.yaml")
DEFAULT_AF2_CONFIG = Path("configs/af2_stage2.yaml")
DEFAULT_SEEDS_CONFIG = Path("configs/seeds.yaml")

PROTEINMPNN_DIR = Path(os.environ.get("PROTEINMPNN_DIR", "/workspace/ProteinMPNN"))
_PMHC_TARGET_CHAINS = ("A", "B", "C")
_BINDER_CHAIN = "D"

# Cycle 03 sub-run A (truncated, 3-chain). After splice_binder_subrun_a the
# spliced complex is A=HLA, B=peptide, C=binder (AF2 FASTA order matching
# LAYOUT_CHAINS["truncated"]); the truncated cleaned target carries HLA on B and
# peptide on C. See PROJECT_STATE 12.j and docs/known_traps.md.
_TRUNCATED_BINDER_CHAIN = "C"
_TRUNCATED_TARGET_CHAINS = ("B", "C")
_VALID_TARGET_LAYOUTS = ("full", "truncated")

HaltStatus = Literal["pass", "halt"]


def _stage1_verdict(stage1: dict[str, Any]) -> str:
    """Resolve Stage 1's pass/fail verdict for the Stage 2 entry gate.

    The cycle-02 merged ``stage1_summary.json`` carries an explicit
    ``halt_rule.verdict``. The cycle-03 sub-run summaries
    (``run_stage1_subrun``) carry ``fraction_geometry_pass`` but no
    ``halt_rule``. Derive the verdict from the *shared* Stage 1 threshold
    (imported, not re-literalled) so Stage 1 and Stage 2 never disagree on
    "did Stage 1 pass".
    """
    halt = stage1.get("halt_rule")
    if isinstance(halt, dict) and "verdict" in halt:
        return str(halt["verdict"])
    frac = stage1.get("fraction_geometry_pass")
    if frac is not None:
        return "PASS" if float(frac) >= STAGE1_HALT_THRESHOLD else "FAIL"
    return "MISSING"


def _resolve_proteinmpnn_python() -> str:
    """Return the SE3nv conda interpreter path for upstream ProteinMPNN.

    Mirrors ``workflow.scripts.run_rfdiffusion._resolve_rfdiffusion_python``:
    pipeline scripts run under uv, but ProteinMPNN's helper scripts live
    in a conda env (torch + the dauparas/ProteinMPNN repo). Overridable
    via ``PROTEINMPNN_PYTHON``; defaults to the SE3nv interpreter that
    bootstrap.sh sets up.
    """
    path = os.environ.get(
        "PROTEINMPNN_PYTHON",
        "/workspace/miniconda3/envs/SE3nv/bin/python",
    )
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"PROTEINMPNN_PYTHON not found: {path}. "
            "Install the SE3nv conda env (see bootstrap.sh) or override "
            "PROTEINMPNN_PYTHON to point at a torch-equipped interpreter."
        )
    return path


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


def _splice_all(
    backbone_pdbs: list[Path],
    cleaned_pmhc: Path,
    spliced_dir: Path,
    target_layout: str = "full",
) -> list[Path]:
    spliced_dir.mkdir(parents=True, exist_ok=True)
    splice_fn = (
        splice_binder.splice_binder_subrun_a
        if target_layout == "truncated"
        else splice_binder.splice_binder_onto_pmhc
    )
    spliced_paths: list[Path] = []
    for backbone in backbone_pdbs:
        out_pdb = spliced_dir / backbone.name
        splice_fn(backbone, cleaned_pmhc, out_pdb)
        spliced_paths.append(out_pdb)
    return spliced_paths


def _extract_target_chain_sequences(
    cleaned_pdb: Path, chains: tuple[str, ...] = _PMHC_TARGET_CHAINS
) -> dict[str, str]:
    """Read the given chains from the cleaned pMHC and return their one-letter
    amino acid sequences, keyed by chain id.

    ``full`` reads A=HC, B=beta2m, C=peptide. ``truncated`` reads B=HLA[1:180]
    and C=peptide from the BAKER groove-only target (no chain A / beta2m).
    """
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(cleaned_pdb.stem, str(cleaned_pdb))
    model = structure[0]
    seqs: dict[str, str] = {}
    for cid in chains:
        if cid not in model:
            raise ValueError(f"{cleaned_pdb}: missing target chain {cid!r}")
        residues = [r for r in model[cid].get_residues() if r.id[0].strip() == ""]
        seqs[cid] = "".join(seq1(r.get_resname()) for r in residues)
    return seqs


def _load_binder_lengths(stage1_designs_jsonl: Path) -> dict[str, int]:
    """Join binder_length by design_id from Stage 1's designs.jsonl
    manifest -- the source of truth for chain D residue counts (never
    recompute from PDBs).
    """
    out: dict[str, int] = {}
    for line in stage1_designs_jsonl.read_text().splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        out[rec["design_id"]] = int(rec["binder_length"])
    return out


def _build_multimer_fasta(
    fasta_path: Path,
    record_id: str,
    target_seqs: dict[str, str],
    binder_seq: str,
    target_layout: str = "full",
) -> None:
    """ColabFold multimer FASTA: a single record with colon-separated chains.
    No trailing colon. Trap #7.

    - ``full`` (cycle 02): order A:B:C:D (HC:beta2m:peptide:binder).
    - ``truncated`` (cycle 03 sub-run A): order HLA:peptide:binder, sourced from
      the truncated cleaned target's chains B (HLA) and C (peptide), so the AF2
      output is A=HLA, B=peptide, C=binder per LAYOUT_CHAINS["truncated"].
    """
    if target_layout == "truncated":
        body = ":".join((target_seqs["B"], target_seqs["C"], binder_seq))
    else:
        body = ":".join((target_seqs["A"], target_seqs["B"], target_seqs["C"], binder_seq))
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    fasta_path.write_text(f">{record_id}\n{body}\n")


def _stage_input_pdbs(spliced_pdbs: list[Path], work_dir: Path) -> Path:
    """Symlink the spliced PDBs into a fresh ``input_pdbs/`` directory.

    ProteinMPNN's ``parse_multiple_chains.py`` greedily reads every file
    under ``--input_path``; staging into a dedicated dir keeps unrelated
    artefacts (Snakemake markers, .ipynb_checkpoints) out of the batch.
    """
    input_dir = work_dir / "input_pdbs"
    if input_dir.exists():
        shutil.rmtree(input_dir)
    input_dir.mkdir(parents=True)
    for pdb in spliced_pdbs:
        link = input_dir / pdb.name
        link.symlink_to(pdb.resolve())
    return input_dir


def _run_real_proteinmpnn(
    spliced_pdbs: list[Path],
    work_dir: Path,
    cfg: dict[str, Any],
    seed: int,
    binder_chain: str = _BINDER_CHAIN,
) -> list[Path]:
    """Drive upstream dauparas/ProteinMPNN via the canonical 3-step
    multi-PDB batch pattern. One model load for the whole batch.

    ``binder_chain`` is the designed chain passed to ``assign_fixed_chains.py
    --chain_list`` (everything else is held fixed): "D" for the cycle-02 full
    layout, "C" for the cycle-03 sub-run A truncated layout (HLA=A, peptide=B,
    binder=C).

    Returns per-design ``.fa`` paths under ``mpnn_out/seqs/`` ordered to
    match ``spliced_pdbs``.
    """
    proteinmpnn_python = _resolve_proteinmpnn_python()
    input_dir = _stage_input_pdbs(spliced_pdbs, work_dir)
    parsed_jsonl = work_dir / "parsed_pdbs.jsonl"
    chains_jsonl = work_dir / "chain_assignments.jsonl"
    mpnn_out = work_dir / "mpnn_out"
    mpnn_out.mkdir(parents=True, exist_ok=True)

    helpers = PROTEINMPNN_DIR / "helper_scripts"
    parse_cmd = [
        proteinmpnn_python,
        str(helpers / "parse_multiple_chains.py"),
        "--input_path",
        str(input_dir),
        "--output_path",
        str(parsed_jsonl),
    ]
    logger.info("ProteinMPNN parse: %s", parse_cmd)
    subprocess.run(parse_cmd, check=True, env={**os.environ})

    assign_cmd = [
        proteinmpnn_python,
        str(helpers / "assign_fixed_chains.py"),
        "--input_path",
        str(parsed_jsonl),
        "--output_path",
        str(chains_jsonl),
        "--chain_list",
        binder_chain,
    ]
    logger.info("ProteinMPNN assign: %s", assign_cmd)
    subprocess.run(assign_cmd, check=True, env={**os.environ})

    run_cmd = [
        proteinmpnn_python,
        str(PROTEINMPNN_DIR / "protein_mpnn_run.py"),
        "--jsonl_path",
        str(parsed_jsonl),
        "--chain_id_jsonl",
        str(chains_jsonl),
        "--out_folder",
        str(mpnn_out),
        "--num_seq_per_target",
        str(cfg["n_seqs_per_backbone"]),
        "--sampling_temp",
        str(cfg["sampling_temperature"]),
        "--seed",
        str(seed),
        "--batch_size",
        "1",
        "--save_score",
        "1",
    ]
    # Cycle 03: penalize RFdiffusion's Ala-heavy compositional prior (Trap
    # #30) via --bias_AA_jsonl when the MPNN config points at a bias file.
    bias_jsonl = cfg.get("bias_AA_jsonl")
    if bias_jsonl:
        run_cmd += ["--bias_AA_jsonl", str(bias_jsonl)]
    logger.info("ProteinMPNN run: %s", run_cmd)
    subprocess.run(run_cmd, check=True, env={**os.environ})

    seqs_dir = mpnn_out / "seqs"
    fastas: list[Path] = []
    for pdb in spliced_pdbs:
        fa = seqs_dir / f"{pdb.stem}.fa"
        if not fa.exists():
            raise FileNotFoundError(f"ProteinMPNN output missing for {pdb.stem}: {fa}")
        fastas.append(fa)
    return fastas


def _run_real_colabfold(
    top_records: list[dict[str, Any]],
    target_seqs: dict[str, str],
    fasta_dir: Path,
    predictions_dir: Path,
    cfg: dict[str, Any],
    target_layout: str = "full",
) -> None:
    """Write one multimer FASTA per ranked record and invoke
    ``colabfold_batch`` for the whole batch. No --use-gpu-relax
    (trap #4).
    """
    fasta_dir.mkdir(parents=True, exist_ok=True)
    predictions_dir.mkdir(parents=True, exist_ok=True)
    ids: list[str] = []
    for rec in top_records:
        record_id = f"{rec['design_id']}_seq{int(rec['seq_id']):02d}"
        fasta_path = fasta_dir / f"{record_id}.fasta"
        _build_multimer_fasta(
            fasta_path, record_id, target_seqs, str(rec["seq"]), target_layout=target_layout
        )
        ids.append(record_id)

    cmd = [
        "colabfold_batch",
        "--num-recycle",
        str(cfg["num_recycles"]),
        "--model-type",
        "alphafold2_multimer_v3",
        "--msa-mode",
        "mmseqs2_uniref",
        "--num-models",
        "5",
        "--rank",
        "multimer",
        str(fasta_dir),
        str(predictions_dir),
    ]
    logger.info("colabfold_batch: %s", cmd)
    subprocess.run(cmd, check=True, env={**os.environ})
    _reshape_flat_outputs(predictions_dir, ids)


def _copy_mock_mpnn_outputs(
    backbones: list[Path], per_design_root: Path, mock_fixtures_dir: Path = MOCK_FIXTURES_DIR
) -> list[Path]:
    src_dir = mock_fixtures_dir / "mpnn_outputs"
    fastas: list[Path] = []
    for backbone in backbones:
        dest_dir = per_design_root / backbone.stem
        dest_dir.mkdir(parents=True, exist_ok=True)
        src = src_dir / f"{backbone.stem}.fa"
        dest = dest_dir / f"{backbone.stem}.fa"
        shutil.copy(src, dest)
        fastas.append(dest)
    return fastas


def _copy_mock_af2_predictions(
    top_records: list[dict[str, Any]],
    predictions_root: Path,
    mock_fixtures_dir: Path = MOCK_FIXTURES_DIR,
) -> None:
    src_root = mock_fixtures_dir / "af2_predictions"
    predictions_root.mkdir(parents=True, exist_ok=True)
    for rec in top_records:
        pred_id = f"{rec['design_id']}_seq{rec['seq_id']}"
        src = src_root / pred_id
        dest = predictions_root / pred_id
        if dest.exists():
            shutil.rmtree(dest)
        shutil.copytree(src, dest)


def _compute_metrics_for_predictions(
    top_records: list[dict[str, Any]],
    predictions_root: Path,
    mock: bool,
    target_layout: str = "full",
) -> list[dict[str, Any]]:
    """Compute iPAE / ipLDDT / BSA for each fan-in prediction.

    Real-mode routes through ``compute_metrics.metrics_for_design``
    (canonical heavy-atom 8 A interface). Mock-mode keeps
    ``metrics_for_control`` (the existing position-slice path) so the
    boundary fixture in ``tests/fixtures/stage2/designs/`` continues to
    yield 4/40 -- see specs/stage2_proteinmpnn_af2.md "Mock mode".
    ``target_layout`` ("full" | "truncated") selects the binder/peptide/MHC
    chain assignment for sub-run A groove-only AF2 outputs (LAYOUT_CHAINS).
    """
    pred_id_fmt = "{design_id}_seq{seq_id}" if mock else "{design_id}_seq{seq_id:02d}"
    metric_records: list[dict[str, Any]] = []
    for rec in top_records:
        pred_id = pred_id_fmt.format(design_id=rec["design_id"], seq_id=int(rec["seq_id"]))
        pred_dir = predictions_root / pred_id
        if mock:
            m = compute_metrics.metrics_for_control(pred_dir, None, target_layout=target_layout)
        else:
            m = compute_metrics.metrics_for_design(pred_dir, target_layout=target_layout)
        record = {
            "design_id": rec["design_id"],
            "seq_id": rec["seq_id"],
            "iPAE": m["ipae"],
            "ipLDDT": m["iplddt"],
            "BSA": m["bsa"],
        }
        # Decomposed iPAE sub-metrics are only produced by the canonical
        # interface path (metrics_for_design); the mock control path omits them.
        if "ppi_pae_int_peptide" in m:
            record["ppi_pae_int_peptide"] = m["ppi_pae_int_peptide"]
            record["ppi_pae_int_mhc"] = m["ppi_pae_int_mhc"]
        metric_records.append(record)
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
    target_layout: str = "full",
    cleaned_pmhc_override: Path | None = None,
) -> int:
    started_at = time.time()
    if target_layout not in _VALID_TARGET_LAYOUTS:
        logger.error(
            "invalid target_layout %r (expected one of %s)", target_layout, _VALID_TARGET_LAYOUTS
        )
        return 1
    binder_chain = _TRUNCATED_BINDER_CHAIN if target_layout == "truncated" else _BINDER_CHAIN

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
    if _stage1_verdict(stage1) != "PASS":
        logger.error("Stage 1 verdict is not PASS at %s; cannot proceed", stage1_summary_path)
        return 1

    target = _load_target(target_manifest)
    if cleaned_pmhc_override is not None:
        cleaned_pmhc = cleaned_pmhc_override
    else:
        cleaned_pmhc = _resolve_cleaned_pdb(target.primary.cleaned_pdb, target_manifest)
    af2_cfg = yaml.safe_load(af2_config.read_text())
    mpnn_cfg = yaml.safe_load(mpnn_config.read_text())
    seeds = aggregate_mpnn_outputs._load_seeds(seeds_config)

    # Locate Stage 1 backbones. In mock mode the fixture ships its own
    # stage1_backbones/ directory next to stage1_summary.json; in real mode
    # backbones live under <stage1_summary_dir>/designs/.
    if mock:
        backbone_dir = stage1_summary_path.parent / "stage1_backbones"
    else:
        backbone_dir = stage1_summary_path.parent / "designs"
    all_backbones = sorted(backbone_dir.glob("design_*.pdb"))
    if not all_backbones:
        logger.error("no stage 1 backbones found under %s", backbone_dir)
        return 1

    # Pre-filter on geometry_pass before MPNN/splice fan-in. Without
    # this, cycle 02 would burn compute on ~170 zero-contact designs
    # and pollute AF2's top-N selection (Trap #29).
    stage1_designs_jsonl = stage1_summary_path.parent / "designs.jsonl"
    binder_lengths: dict[str, int] = {}
    if stage1_designs_jsonl.exists():
        all_records = [
            json.loads(line)
            for line in stage1_designs_jsonl.read_text().splitlines()
            if line.strip()
        ]
        pass_ids = {r["design_id"] for r in all_records if r.get("geometry_pass")}
        backbones = [b for b in all_backbones if b.stem in pass_ids]
        logger.info(
            "Stage 2 input filtered: %d/%d Stage 1 designs pass geometry_pass=True",
            len(backbones),
            len(all_backbones),
        )
        if not backbones:
            logger.error(
                "no geometry-pass designs available for Stage 2; halt "
                "(checked %d Stage 1 records)",
                len(all_records),
            )
            return 1
        if not mock:
            binder_lengths = _load_binder_lengths(stage1_designs_jsonl)
    elif mock:
        # Mock fixture variant without a designs.jsonl: trust the glob.
        backbones = all_backbones
    else:
        logger.error("stage1 designs manifest missing: %s", stage1_designs_jsonl)
        return 1

    out_dir.mkdir(parents=True, exist_ok=True)
    spliced_dir = out_dir / "proteinmpnn_designs" / "spliced"
    spliced_pdbs = _splice_all(backbones, cleaned_pmhc, spliced_dir, target_layout=target_layout)
    (spliced_dir / ".done").touch()

    # Mock MPNN/AF2 fixtures are co-located with the Stage 1 backbones the
    # mock run consumed, so the truncated fixture set (3-chain predictions) is
    # used for the truncated path rather than the 4-chain full-layout default.
    mock_fixtures_dir = stage1_summary_path.parent
    per_design_root = out_dir / "proteinmpnn_designs" / "per_design"
    if mock:
        per_design_fastas = _copy_mock_mpnn_outputs(backbones, per_design_root, mock_fixtures_dir)
    else:
        # Cycle 02 seed contract: a single MPNN seed governs the whole
        # batch (one-call invocation; per-design / per-sample seeding is
        # deferred to cycle 03). We assert the first reserved seed for
        # this cycle exists, matching aggregator behaviour.
        batch_seed = aggregate_mpnn_outputs._compute_mpnn_seed(cycle, 0, 0, seeds)
        mpnn_work_dir = out_dir / "proteinmpnn_designs" / "mpnn_work"
        mpnn_work_dir.mkdir(parents=True, exist_ok=True)
        per_design_fastas = _run_real_proteinmpnn(
            spliced_pdbs, mpnn_work_dir, mpnn_cfg, batch_seed, binder_chain=binder_chain
        )

    design_records: list[dict[str, str | int | None]] = [
        {
            "design_id": backbone.stem,
            "backbone_pdb": str(backbone),
            "spliced_pdb": str(spliced),
            "binder_length": binder_lengths.get(backbone.stem) if not mock else None,
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

    # Real-mode rank: ProteinMPNN's global_score (lower=better) per spec.
    # Mock fixtures only emit ``score`` so fall back to mpnn_nll there.
    def _rank_key(rec: dict[str, Any]) -> float:
        if "mpnn_global_score" in rec:
            return float(rec["mpnn_global_score"])
        return float(rec["mpnn_nll"])

    records.sort(key=_rank_key)
    n_total = len(records)
    fan_in_top_n = min(n_total, 40) if mock else min(int(af2_cfg["fan_in_top_n"]), n_total)
    top_records = records[:fan_in_top_n]

    predictions_root = out_dir / "af2_designs" / "predictions"
    if mock:
        _copy_mock_af2_predictions(top_records, predictions_root, mock_fixtures_dir)
    else:
        target_chains = (
            _TRUNCATED_TARGET_CHAINS if target_layout == "truncated" else _PMHC_TARGET_CHAINS
        )
        target_seqs = _extract_target_chain_sequences(cleaned_pmhc, chains=target_chains)
        fasta_dir = out_dir / "af2_designs" / "fastas"
        _run_real_colabfold(
            top_records,
            target_seqs,
            fasta_dir,
            predictions_root,
            af2_cfg,
            target_layout=target_layout,
        )

    metric_records = _compute_metrics_for_predictions(
        top_records, predictions_root, mock=mock, target_layout=target_layout
    )
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
        "--target-layout",
        choices=_VALID_TARGET_LAYOUTS,
        default="full",
        help="full = cycle-02 4-chain (binder D); truncated = cycle-03 sub-run A 3-chain (binder C).",
    )
    parser.add_argument(
        "--cleaned-pmhc",
        type=Path,
        default=None,
        help="Override the cleaned target PDB for splice/FASTA (truncated path supplies "
        "the BAKER groove-only target instead of the manifest's full cleaned_pdb).",
    )
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
        target_layout=args.target_layout,
        cleaned_pmhc_override=args.cleaned_pmhc,
    )


if __name__ == "__main__":
    sys.exit(main())
