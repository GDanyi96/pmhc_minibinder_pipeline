"""Tests for the BAKER-truncated controls path + target-layout metrics adapter.

The truncated controls fold the binder against the BAKER groove-only fragment
(HLA[1:180] : peptide : binder, beta2m dropped), yielding a 3-chain AF2 output
A=HLA, B=peptide, C=binder. compute_metrics must decompose iPAE correctly for
both the full (4-chain) and truncated (3-chain) layouts.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np

import scripts.run_controls as run_controls
from workflow.scripts import compute_metrics

REPO_ROOT = Path(__file__).resolve().parent.parent


# --- target layout chain assignments -----------------------------------------


def test_full_layout_matches_legacy_defaults() -> None:
    full = compute_metrics.LAYOUT_CHAINS["full"]
    assert full["binder"] == "D"
    assert full["peptide"] == ("C",)
    assert full["mhc"] == ("A", "B")
    assert full["combined"] == ("A", "B", "C")


def test_truncated_layout_assignment() -> None:
    trunc = compute_metrics.LAYOUT_CHAINS["truncated"]
    assert trunc["binder"] == "C"
    assert trunc["peptide"] == ("B",)
    assert trunc["mhc"] == ("A",)
    assert trunc["combined"] == ("A", "B")


# --- FASTA composition --------------------------------------------------------


def test_truncated_fasta_drops_beta2m() -> None:
    controls = run_controls.load_manifest(run_controls.DEFAULT_MANIFEST)
    p1 = next(c for c in controls if c.id == "P1")
    _, full_body = run_controls.build_fasta_entry(p1, run_controls.DEFAULT_WT1_PDB, mock=True)
    _, trunc_body = run_controls.build_fasta_entry(
        p1, run_controls.DEFAULT_TRUNCATED_PDB, mock=True, target_layout="truncated"
    )
    assert full_body.count(":") == 3  # hc:b2m:peptide:binder
    assert trunc_body.count(":") == 2  # hla:peptide:binder
    hla, peptide, binder = trunc_body.split(":")
    assert len(hla) == 180  # HLA groove only, no beta2m
    assert peptide == p1.target_peptide
    assert binder == p1.binder_sequence


# --- slice-path metrics (controls) -------------------------------------------


def test_compute_ipae_truncated_binder_chain() -> None:
    # 2 A (HLA), 1 B (peptide), 2 C (binder) -> N=5
    slices = compute_metrics.chain_slices({"A": 2, "B": 1, "C": 2})
    pae = [[1.0] * 5 for _ in range(5)]
    # binder C rows/cols (indices 3,4) vs A (0,1) and B (2): set distinct value
    for d in (3, 4):
        for j in (0, 1, 2):
            pae[d][j] = pae[j][d] = 5.0
    ipae = compute_metrics.compute_ipae(pae, slices, binder="C", mhc_chains=("A", "B"))
    assert ipae == 5.0
    plddt = [10.0, 10.0, 10.0, 80.0, 90.0]
    assert compute_metrics.compute_iplddt(plddt, slices, binder="C") == 85.0


# --- interface-path decomposition (designs) ----------------------------------
# 3-chain truncated layout: A=HLA(2), B=peptide(1), C=binder(2).
_TRUNC_CHAIN_IDS = np.asarray(["A", "A", "B", "C", "C"], dtype="<U1")


def _trunc_distances() -> np.ndarray:
    d = np.full((5, 5), 50.0)
    np.fill_diagonal(d, 0.0)
    d[3, 0] = d[0, 3] = 5.0  # binder C[3] <-> HLA A[0]
    d[4, 2] = d[2, 4] = 6.0  # binder C[4] <-> peptide B[2]
    return d


def _trunc_pae() -> np.ndarray:
    pae = np.full((5, 5), 30.0)
    pae[3, 0] = pae[0, 3] = 2.0
    pae[4, 2] = pae[2, 4] = 4.0
    return pae


def test_truncated_decomposition_uses_correct_chains() -> None:
    chains = compute_metrics.LAYOUT_CHAINS["truncated"]
    pae, dist = _trunc_pae(), _trunc_distances()
    mhc = compute_metrics.compute_ipae_interface(
        pae, _TRUNC_CHAIN_IDS, dist, binder=chains["binder"], target=chains["mhc"]
    )
    peptide = compute_metrics.compute_ipae_interface(
        pae, _TRUNC_CHAIN_IDS, dist, binder=chains["binder"], target=chains["peptide"]
    )
    assert mhc == 2.0  # C <-> A only
    assert peptide == 4.0  # C <-> B only


def test_full_layout_default_misses_3chain_interface() -> None:
    """The full default (binder='D') finds no interface on a 3-chain truncated
    prediction -- the precise failure the target-layout adapter prevents."""
    pae, dist = _trunc_pae(), _trunc_distances()
    full = compute_metrics.compute_ipae_interface(pae, _TRUNC_CHAIN_IDS, dist)  # binder="D"
    assert full == float("inf")
    chains = compute_metrics.LAYOUT_CHAINS["truncated"]
    trunc = compute_metrics.compute_ipae_interface(
        pae, _TRUNC_CHAIN_IDS, dist, binder=chains["binder"], target=chains["combined"]
    )
    assert np.isfinite(trunc)


# --- run() integration: dirs + filename --------------------------------------


def test_run_truncated_mock_writes_baseline(tmp_path: Path) -> None:
    # Uncalibrated mock ColabFold fixtures trip the P2 calibration halt gate
    # (rc=1) -- expected; metrics.json is still written before the halt return.
    cycle = "97"
    results = REPO_ROOT / f"results/cycle_{cycle}"
    if results.exists():
        shutil.rmtree(results)
    try:
        run_controls.main(["--mock", "--target", "truncated", "--cycle", cycle])
        baseline = results / "controls_truncated_baseline"
        metrics = baseline / "metrics_truncated_baseline.json"
        assert metrics.exists()
        doc = json.loads(metrics.read_text())
        assert doc["target"] == "truncated"
        fasta = baseline / "colabfold" / "_controls_fasta" / "P1.fasta"
        assert fasta.read_text().strip().splitlines()[1].count(":") == 2
        # The canonical full-target controls dir must be untouched.
        assert not (results / "stage2" / "metrics.json").exists()
    finally:
        if results.exists():
            shutil.rmtree(results)
