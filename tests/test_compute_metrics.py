"""Unit tests for workflow.scripts.compute_metrics."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from workflow.scripts import compute_metrics


def _write_scores(path: Path, pae_key: str, pae_value: list[list[float]]) -> None:
    path.write_text(json.dumps({pae_key: pae_value, "plddt": [90.0]}))


def test_load_pae_new_key(tmp_path: Path) -> None:
    scores = tmp_path / "scores.json"
    _write_scores(scores, "pae", [[0.0, 1.0], [1.0, 0.0]])
    assert compute_metrics.load_pae(scores) == [[0.0, 1.0], [1.0, 0.0]]


def test_load_pae_legacy_key(tmp_path: Path) -> None:
    scores = tmp_path / "scores.json"
    _write_scores(scores, "predicted_aligned_error", [[2.0, 3.0], [3.0, 2.0]])
    assert compute_metrics.load_pae(scores) == [[2.0, 3.0], [3.0, 2.0]]


def test_load_pae_missing_raises(tmp_path: Path) -> None:
    scores = tmp_path / "scores.json"
    scores.write_text(json.dumps({"plddt": [90.0]}))
    with pytest.raises(KeyError, match="missing both"):
        compute_metrics.load_pae(scores)


# ---------- Canonical interface-mask metrics ----------


def _toy_chain_ids() -> np.ndarray:
    # 2 A, 1 B, 1 C, 2 D = 6 residues
    return np.asarray(["A", "A", "B", "C", "D", "D"], dtype="<U1")


def _toy_distances_with_interface() -> np.ndarray:
    """D[0] is within 8 A of A[0] and C[0]; D[1] is within 8 A of A[1]
    only. All other D<->ABC pairs are > 8 A. Within-chain entries are
    irrelevant to the interface mask but set to 1.0 for realism.
    """
    d = np.full((6, 6), 50.0)
    np.fill_diagonal(d, 0.0)
    # Within-chain
    d[0, 1] = d[1, 0] = 1.0
    d[4, 5] = d[5, 4] = 1.0
    # Interface pairs (D rows 4, 5)
    d[4, 0] = d[0, 4] = 5.0
    d[4, 3] = d[3, 4] = 7.5
    d[5, 1] = d[1, 5] = 6.0
    # All other D-vs-ABC kept at 50.0
    return d


def test_compute_ipae_interface_known_fixture() -> None:
    chain_ids = _toy_chain_ids()
    distances = _toy_distances_with_interface()
    # PAE matrix: asymmetric on the interface pairs so the min(pae_ij,
    # pae_ji) reduction is exercised. Non-interface entries set to a
    # large value to surface bugs if the mask is wrong.
    pae = np.full((6, 6), 99.0)
    pae[4, 0], pae[0, 4] = 4.0, 6.0  # min = 4.0
    pae[4, 3], pae[3, 4] = 8.0, 7.0  # min = 7.0
    pae[5, 1], pae[1, 5] = 10.0, 9.0  # min = 9.0
    ipae = compute_metrics.compute_ipae_interface(pae, chain_ids, distances)
    # Three interface pairs in the D-rows x ABC-cols mask: (4,0), (4,3), (5,1).
    # Mean of 4.0, 7.0, 9.0 = 20/3.
    assert ipae == pytest.approx(20.0 / 3.0)


def test_compute_iplddt_interface_known_fixture() -> None:
    chain_ids = _toy_chain_ids()
    distances = _toy_distances_with_interface()
    # Residues touched by the interface mask: D rows 4, 5 (both have
    # interface pairs); ABC cols 0, 1, 3 (A[0], A[1], C[0]). B[0]=index 2
    # is not in the interface. Expected = mean over indices {0, 1, 3, 4, 5}.
    plddt = np.asarray([90.0, 92.0, 50.0, 88.0, 70.0, 72.0])
    iplddt = compute_metrics.compute_iplddt_interface(plddt, chain_ids, distances)
    expected = float(np.mean([90.0, 92.0, 88.0, 70.0, 72.0]))
    assert iplddt == pytest.approx(expected)


def test_no_interface_yields_inf_ipae_and_zero_iplddt() -> None:
    chain_ids = _toy_chain_ids()
    distances = np.full((6, 6), 50.0)
    np.fill_diagonal(distances, 0.0)
    pae = np.full((6, 6), 4.0)
    plddt = np.full(6, 90.0)
    assert compute_metrics.compute_ipae_interface(pae, chain_ids, distances) == float("inf")
    assert compute_metrics.compute_iplddt_interface(plddt, chain_ids, distances) == 0.0


def test_compute_heavy_atom_distances_excludes_hydrogens(tmp_path: Path) -> None:
    """A minimal 2-residue PDB with one heavy + one H atom per residue:
    residue-residue min distance must be the heavy-heavy distance, not
    the closer H-H distance.
    """
    pdb = tmp_path / "mini.pdb"
    # Residue 1: CA at (0,0,0), H at (10,0,0). Residue 2: CA at (4,0,0),
    # H at (10.5,0,0). Heavy-heavy min = 4.0. H-H min = 0.5 (must be
    # ignored by the distance routine).
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C\n"
        "ATOM      2  H   ALA A   1      10.000   0.000   0.000  1.00 90.00           H\n"
        "ATOM      3  CA  ALA A   2       4.000   0.000   0.000  1.00 90.00           C\n"
        "ATOM      4  H   ALA A   2      10.500   0.000   0.000  1.00 90.00           H\n"
        "END\n"
    )
    distances = compute_metrics.compute_heavy_atom_distances(pdb)
    assert distances.shape == (2, 2)
    assert distances[0, 1] == pytest.approx(4.0)
    assert distances[0, 0] == 0.0


def test_load_plddt_from_pdb_reads_ca_bfactor(tmp_path: Path) -> None:
    pdb = tmp_path / "plddt.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 87.50           C\n"
        "ATOM      2  CA  ALA A   2       3.800   0.000   0.000  1.00 93.20           C\n"
        "END\n"
    )
    plddt = compute_metrics.load_plddt_from_pdb(pdb)
    assert plddt.tolist() == pytest.approx([87.50, 93.20])


def test_chain_ids_per_residue_matches_pdb_order(tmp_path: Path) -> None:
    pdb = tmp_path / "chains.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C\n"
        "ATOM      2  CA  ALA A   2       3.800   0.000   0.000  1.00 90.00           C\n"
        "ATOM      3  CA  ALA B   1      20.000   0.000   0.000  1.00 90.00           C\n"
        "ATOM      4  CA  ALA D   1       0.000   0.000   5.000  1.00 90.00           C\n"
        "END\n"
    )
    ids = compute_metrics.chain_ids_per_residue(pdb)
    assert ids.tolist() == ["A", "A", "B", "D"]
