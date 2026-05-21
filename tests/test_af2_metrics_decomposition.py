"""Tests for the cycle-03 decomposed iPAE sub-metrics in compute_metrics.

BAKER_LAB_2025's ppi_pae_int is our combined iPAE (binder D <-> A+B+C). Cycle
03 splits it into peptide-facing (D <-> C) and MHC-facing (D <-> A+B) parts.
"""

from __future__ import annotations

import numpy as np

from workflow.scripts import compute_metrics

# 2 A, 1 B, 1 C, 2 D = 6 residues
CHAIN_IDS = np.asarray(["A", "A", "B", "C", "D", "D"], dtype="<U1")


def _distances() -> np.ndarray:
    d = np.full((6, 6), 50.0)
    np.fill_diagonal(d, 0.0)
    # D[4] within 8 A of A[0] (MHC) and C[3] (peptide); D[5] near A[1] only.
    d[4, 0] = d[0, 4] = 5.0
    d[4, 3] = d[3, 4] = 7.0
    d[5, 1] = d[1, 5] = 6.0
    return d


def _pae() -> np.ndarray:
    pae = np.full((6, 6), 30.0)
    pae[4, 0] = pae[0, 4] = 2.0
    pae[4, 3] = pae[3, 4] = 4.0
    pae[5, 1] = pae[1, 5] = 6.0
    return pae


def test_peptide_submetric_uses_only_chain_c() -> None:
    pae, dist = _pae(), _distances()
    peptide = compute_metrics.compute_ipae_interface(pae, CHAIN_IDS, dist, target=("C",))
    # Only the D[4]<->C[3] pair (pae 4.0) qualifies.
    assert peptide == 4.0


def test_mhc_submetric_uses_a_and_b() -> None:
    pae, dist = _pae(), _distances()
    mhc = compute_metrics.compute_ipae_interface(pae, CHAIN_IDS, dist, target=("A", "B"))
    # D[4]<->A[0] (2.0) and D[5]<->A[1] (6.0) -> mean 4.0.
    assert mhc == 4.0


def test_combined_equals_default_target() -> None:
    pae, dist = _pae(), _distances()
    combined = compute_metrics.compute_ipae_interface(pae, CHAIN_IDS, dist)
    default = compute_metrics.compute_ipae_interface(pae, CHAIN_IDS, dist, target=("A", "B", "C"))
    assert combined == default


def test_no_peptide_interface_returns_inf() -> None:
    pae = _pae()
    dist = _distances()
    dist[4, 3] = dist[3, 4] = 50.0  # break the only peptide contact
    peptide = compute_metrics.compute_ipae_interface(pae, CHAIN_IDS, dist, target=("C",))
    assert peptide == float("inf")
