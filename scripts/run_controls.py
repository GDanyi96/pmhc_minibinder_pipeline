# mypy: disable-error-code="no-untyped-call,attr-defined,no-any-return,misc"
"""Stage 2 — Part A controls orchestrator.

Drives the controls validation flow:
1. Load manifest + thresholds (pydantic).
2. Assert P1 sequence integrity (I12, E40 — BAKER_LAB_2025 fingerprint).
3. Populate N1/N2 via generate_negatives if needed.
4. Write per-control FASTA files; run ColabFold (or mock).
5. Compute iPAE/ipLDDT/BSA per control.
6. Enforce halt gate (P1 absolute, P2 calibration, P1 rank, P3 specificity).
7. Emit results/cycle_NN/stage2/metrics.json.

Halt rules and thresholds live in configs/thresholds.yaml under stage2.controls.
"""

from __future__ import annotations

import os

# RunPod base sets LD_LIBRARY_PATH=/usr/local/cuda/lib64 (CUDA 13), which
# overrides JAX's bundled CUDA 12 libs and SIGSEGVs cuDNN. The colabfold
# subprocess we spawn inherits this process's env, so popping here
# propagates the unset to the JAX child. See bootstrap.sh.
os.environ.pop("LD_LIBRARY_PATH", None)

import argparse
import datetime as dt
import json
import logging
import sys
from collections.abc import Iterable
from pathlib import Path
from typing import Any, Literal

import yaml
from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
from pydantic import BaseModel, ConfigDict

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts import generate_negatives
from workflow.scripts import compute_metrics, run_colabfold

logger = logging.getLogger("run_controls")
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

THREE_TO_ONE: dict[str, str] = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

DEFAULT_MANIFEST = Path("data/controls/controls_manifest.yaml")
DEFAULT_THRESHOLDS = Path("configs/thresholds.yaml")
DEFAULT_CYCLE = "01"
DEFAULT_WT1_PDB = Path("data/targets/3hpj_clean.pdb")


class ControlEntry(BaseModel):
    model_config = ConfigDict(extra="allow")
    id: str
    label: str
    binder_sequence: str | None = None
    binder_length: int
    target_peptide: str
    target_hla: str
    target_pdb: str | None = None
    backbone_pdb: str | None = None
    halt_on_fail: bool = False
    expected_ipae_max: float | None = None
    expected_ipae: float | None = None
    expected_ipae_tolerance: float | None = None
    expected_ipae_min_delta: float | None = None
    expected_ipae_min: float | None = None


class ControlsThresholds(BaseModel):
    p1_ipae_max: float
    p2_ipae_tolerance: float
    pos_neg_ipae_gap_min: float
    pos_neg_iplddt_gap_min: float


class Stage2Thresholds(BaseModel):
    controls: ControlsThresholds


class Thresholds(BaseModel):
    stage2: Stage2Thresholds


def load_manifest(path: Path) -> list[ControlEntry]:
    with path.open() as fh:
        raw = yaml.safe_load(fh)
    return [ControlEntry.model_validate(e) for e in raw["controls"]]


def load_thresholds(path: Path) -> ControlsThresholds:
    with path.open() as fh:
        raw = yaml.safe_load(fh)
    return Thresholds.model_validate(raw).stage2.controls


def assert_p1_integrity(p1: ControlEntry) -> None:
    """BAKER_LAB_2025 structural fingerprint anchors: I12 H-bond to peptide Y8,
    E40 bidentate to peptide R1. A silent edit to wt1-5 must not propagate.
    """
    seq = p1.binder_sequence
    if not seq:
        raise SystemExit("P1 sequence integrity check failed: binder_sequence is empty")
    if len(seq) < 40 or seq[11] != "I" or seq[39] != "E":
        raise SystemExit(
            "P1 sequence integrity check failed at I12/E40 — manifest corrupted "
            f"(got seq[11]={seq[11] if len(seq) > 11 else '?'!r}, "
            f"seq[39]={seq[39] if len(seq) > 39 else '?'!r})"
        )


def chain_sequence(pdb_path: Path, chain_id: str) -> str:
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = structure[0]
    if chain_id not in {c.id for c in model.get_chains()}:
        raise ValueError(f"{pdb_path}: chain {chain_id} absent")
    chain = model[chain_id]
    residues = [r for r in chain.get_residues() if r.id[0].strip() == ""]
    return "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in residues)


def build_fasta_entry(
    control: ControlEntry,
    wt1_pdb: Path,
    mock: bool = False,
) -> tuple[str, str]:
    backbone = control.backbone_pdb
    if mock:
        # Mock mode: real HC/beta2m chain extraction is unnecessary — the mock
        # ColabFold wrapper consumes only the FASTA header. Emit length-correct
        # placeholder sequences so real and mock FASTAs have the same shape.
        hc, b2m = "G" * 275, "G" * 100
    else:
        # P3 (backbone_pdb=null) uses HLA-A*02:01 HC + human beta2m from the WT1
        # cleaned PDB (same allele as P3) per HOUSEHOLDER_GARCIA_2025.
        source_pdb = wt1_pdb if backbone is None else Path(backbone)
        hc = chain_sequence(source_pdb, "A")
        b2m = chain_sequence(source_pdb, "B")
    peptide = control.target_peptide
    binder = control.binder_sequence
    if not binder:
        raise ValueError(f"control {control.id}: binder_sequence missing")
    body = f"{hc}:{b2m}:{peptide}:{binder}"
    return control.id, body


def write_fastas(
    controls: Iterable[ControlEntry],
    out_dir: Path,
    wt1_pdb: Path,
    mock: bool = False,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for control in controls:
        ident, body = build_fasta_entry(control, wt1_pdb, mock=mock)
        fasta = out_dir / f"{ident}.fasta"
        fasta.write_text(f">{ident}\n{body}\n")


HaltStatus = Literal["pass", "halt"]


def enforce_halt_gate(
    metrics: dict[str, dict[str, Any]],
    thresholds: ControlsThresholds,
) -> tuple[HaltStatus, str]:
    """Biology-correct controls halt gate.

    Replaces the cycle-1 "P1 rank <= 2 of 5" rule, which fired falsely on
    real data: P2 has more interface (longer published binder) and P3 ~ P1
    due to the AF2-multimer specificity blind spot (HOUSEHOLDER_GARCIA_2025,
    MARES_IOANNIDIS_2025). Specificity is Stage 3 cross-pan's job.

    Four checks: P1 cognate binding, P2 calibration vs Jenkins published
    value, positives <-> negatives iPAE gap, positives <-> negatives
    ipLDDT gap.
    """
    p1, p2, p3 = (metrics[k]["ipae"] for k in ("P1", "P2", "P3"))
    n1, n2 = (metrics[k]["ipae"] for k in ("N1", "N2"))
    p1l, p2l, p3l = (metrics[k]["iplddt"] for k in ("P1", "P2", "P3"))
    n1l, n2l = (metrics[k]["iplddt"] for k in ("N1", "N2"))

    if p1 > thresholds.p1_ipae_max:
        return "halt", f"HALT: P1 iPAE {p1:.2f} > {thresholds.p1_ipae_max:.1f} Ang threshold"

    if abs(p2 - 6.5) > thresholds.p2_ipae_tolerance:
        return "halt", (
            f"HALT: ColabFold calibration failure - P2 iPAE {p2:.2f} expected "
            f"6.5 +/- {thresholds.p2_ipae_tolerance:.1f} Ang. "
            "Check chain ordering and num_recycles."
        )

    ipae_gap = min(n1, n2) - max(p1, p2, p3)
    if ipae_gap < thresholds.pos_neg_ipae_gap_min:
        return "halt", (
            f"HALT: positives<->negatives iPAE gap {ipae_gap:.2f} Ang "
            f"< {thresholds.pos_neg_ipae_gap_min:.1f} Ang "
            f"(pos max {max(p1, p2, p3):.2f}, neg min {min(n1, n2):.2f})"
        )

    iplddt_gap = min(p1l, p2l, p3l) - max(n1l, n2l)
    if iplddt_gap < thresholds.pos_neg_iplddt_gap_min:
        return "halt", (
            f"HALT: positives<->negatives ipLDDT gap {iplddt_gap:.2f} "
            f"< {thresholds.pos_neg_iplddt_gap_min:.1f} "
            f"(pos min {min(p1l, p2l, p3l):.2f}, neg max {max(n1l, n2l):.2f})"
        )

    logger.info(
        "Stage 2 does not assess specificity - P3 (binder vs non-cognate pMHC) "
        "iPAE %.2f Ang, P1 (binder vs cognate) iPAE %.2f Ang. AF2-multimer "
        "specificity blind spot per HOUSEHOLDER_GARCIA_2025, MARES_IOANNIDIS_2025. "
        "Specificity is assessed in Stage 3 cross-panning, not here.",
        p3,
        p1,
    )
    return "pass", ""


def _passes_per_control(control: ControlEntry, ipae: float) -> bool:
    if control.expected_ipae_max is not None:
        return ipae <= control.expected_ipae_max
    if control.expected_ipae is not None and control.expected_ipae_tolerance is not None:
        return abs(ipae - control.expected_ipae) <= control.expected_ipae_tolerance
    if control.expected_ipae_min is not None:
        return ipae >= control.expected_ipae_min
    return True


def run(
    manifest_path: Path,
    thresholds_path: Path,
    cycle: str,
    mock: bool,
    wt1_pdb: Path,
    skip_colabfold: bool = False,
    metrics_only: bool = False,
) -> int:
    controls = load_manifest(manifest_path)
    thresholds = load_thresholds(thresholds_path)

    by_id = {c.id: c for c in controls}
    assert_p1_integrity(by_id["P1"])

    if not metrics_only and any(by_id[i].binder_sequence is None for i in ("N1", "N2")):
        logger.info("Populating N1/N2 via generate_negatives (seed=42)")
        generate_negatives.main(["--config", str(manifest_path)])
        controls = load_manifest(manifest_path)
        by_id = {c.id: c for c in controls}

    results_dir = Path(f"results/cycle_{cycle}/stage2")
    fasta_dir = results_dir / "colabfold" / "_controls_fasta"
    colabfold_out = results_dir / "colabfold" / "controls"
    if not metrics_only:
        write_fastas(controls, fasta_dir, wt1_pdb, mock=mock)

    if not (skip_colabfold or metrics_only):
        cf_argv = ["--fasta-dir", str(fasta_dir), "--out-dir", str(colabfold_out)]
        if mock:
            cf_argv.append("--mock")
        run_colabfold.main(cf_argv)
    else:
        logger.info(
            "Skipping ColabFold (%s); reusing existing outputs under %s",
            "--metrics-only" if metrics_only else "--skip-colabfold",
            colabfold_out,
        )

    per_control: dict[str, dict[str, Any]] = {}
    for control in controls:
        directory = colabfold_out / control.id
        row = compute_metrics.metrics_for_control(directory, None)
        per_control[control.id] = row

    status, halt_message = enforce_halt_gate(per_control, thresholds)

    for control in controls:
        row = per_control[control.id]
        row["pass"] = _passes_per_control(control, row["ipae"])

    metrics_doc = {
        "cycle": cycle,
        "stage": "2",
        "timestamp": dt.datetime.now(dt.UTC).isoformat(),
        "controls": per_control,
        "halt_gate": {"passed": status == "pass", "message": halt_message},
        "designs": [],
    }

    metrics_path = results_dir / "metrics.json"
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    metrics_path.write_text(json.dumps(metrics_doc, indent=2) + "\n")
    logger.info("Wrote %s", metrics_path)

    if status == "halt":
        logger.error(halt_message)
        return 1
    logger.info(
        "Controls halt gate passed (P1 cognate, P2 calibration, positives<->negatives gaps)."
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--thresholds", type=Path, default=DEFAULT_THRESHOLDS)
    parser.add_argument("--cycle", default=DEFAULT_CYCLE)
    parser.add_argument("--wt1-pdb", type=Path, default=DEFAULT_WT1_PDB)
    parser.add_argument("--mock", action="store_true")
    parser.add_argument(
        "--skip-colabfold",
        action="store_true",
        help="Assume ColabFold outputs already exist; skip the subprocess call.",
    )
    parser.add_argument(
        "--metrics-only",
        action="store_true",
        help="Skip negatives gen + FASTA writing + ColabFold; recompute metrics+halt only.",
    )
    args = parser.parse_args(argv)
    return run(
        args.config,
        args.thresholds,
        args.cycle,
        args.mock,
        args.wt1_pdb,
        skip_colabfold=args.skip_colabfold,
        metrics_only=args.metrics_only,
    )


if __name__ == "__main__":
    sys.exit(main())
