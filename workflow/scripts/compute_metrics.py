# mypy: disable-error-code="no-untyped-call,attr-defined,no-any-return,misc"
"""Stage 2c — iPAE, ipLDDT, and BSA computation from ColabFold outputs.

Two metric flavors live here:

* Position-slice iPAE/ipLDDT (``compute_ipae``/``compute_iplddt``/
  ``metrics_for_control``) — the cycle-1 controls' calibrated path.
  iPAE is the symmetric mean PAE over the full D <-> A+B+C submatrix;
  ipLDDT is the mean pLDDT over chain D residues. Fast, no PDB
  geometry, no interface threshold.

* Canonical heavy-atom interface iPAE/ipLDDT
  (``compute_ipae_interface``/``compute_iplddt_interface``/
  ``metrics_for_design``) — HADRUP_JENKINS_2025 Fig 1B operational
  definition: iPAE = mean min(pae_ij, pae_ji) over residue pairs where
  one residue is on chain D, the other on chains A/B/C, and the
  heavy-atom min distance is <= 8 A. ipLDDT = mean pLDDT over residues
  that touch the interface from either side. The new Stage 2 designs
  pipeline (``scripts/run_stage2.py``) routes real-mode predictions
  through this canonical path; the controls pipeline keeps the slice
  path to preserve cycle-1 calibration.

BSA: SASA(pMHC) + SASA(binder) - SASA(complex) via freesasa default
Lee-Richards. Splits the complex PDB into A+B+C and D-only temporary PDBs
written via Biopython PDBIO so freesasa.Structure.from_file consumes only
ATOM records (pitfall #8).

Chain order A=HC, B=beta2m, C=peptide, D=binder. Boundaries are computed
from the actual PDB at runtime, never hard-coded (pitfall #4). The JSON
key `iptm` is never used as an iPAE proxy (pitfall #2).
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from collections.abc import Sequence
from glob import glob
from pathlib import Path
from typing import Any

import freesasa
import numpy as np
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure

CHAIN_ORDER = ("A", "B", "C", "D")
DEFAULT_SCORES_GLOB = "*_scores_rank_001*.json"
DEFAULT_PDB_GLOB = "*_unrelaxed_rank_001*.pdb"


def load_pae(scores_json: Path) -> list[list[float]]:
    payload = json.loads(scores_json.read_text())
    # ColabFold >= 1.6.0 renamed the key to "pae"; older versions used
    # "predicted_aligned_error". Prefer the new name; fall back for compat
    # with legacy fixtures and pre-1.6 ColabFold installs.
    pae = payload.get("pae")
    if pae is None:
        pae = payload.get("predicted_aligned_error")
    if pae is None:
        raise KeyError(f"{scores_json}: missing both 'pae' and 'predicted_aligned_error'")
    return pae


def load_plddt(scores_json: Path) -> list[float]:
    payload = json.loads(scores_json.read_text())
    plddt = payload.get("plddt")
    if plddt is None:
        raise KeyError(f"{scores_json}: missing 'plddt'")
    return plddt


def chain_lengths(pdb_path: Path) -> dict[str, int]:
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = structure[0]
    lengths: dict[str, int] = {}
    for chain in model.get_chains():
        if chain.id in CHAIN_ORDER:
            lengths[chain.id] = sum(1 for r in chain.get_residues() if r.id[0].strip() == "")
    return lengths


def chain_slices(lengths: dict[str, int]) -> dict[str, slice]:
    slices: dict[str, slice] = {}
    cursor = 0
    for chain_id in CHAIN_ORDER:
        size = lengths.get(chain_id, 0)
        slices[chain_id] = slice(cursor, cursor + size)
        cursor += size
    return slices


def compute_ipae(pae: list[list[float]], slices: dict[str, slice]) -> float:
    d_range = range(slices["D"].start, slices["D"].stop)
    mhc_ranges: list[range] = [range(slices[c].start, slices[c].stop) for c in ("A", "B", "C")]
    total = 0.0
    count = 0
    for d in d_range:
        for r in mhc_ranges:
            for j in r:
                total += pae[d][j] + pae[j][d]
                count += 2
    if count == 0:
        raise ValueError("compute_ipae: empty cross-chain submatrix (check chain ordering)")
    return total / count


def compute_iplddt(plddt: list[float], slices: dict[str, slice]) -> float:
    binder = plddt[slices["D"]]
    if not binder:
        raise ValueError("compute_iplddt: chain D is empty (check chain ordering)")
    return sum(binder) / len(binder)


class _DropChain(Select):
    def __init__(self, drop: str) -> None:
        self._drop = drop

    def accept_chain(self, chain: Chain) -> bool:
        return chain.id != self._drop


class _KeepChain(Select):
    def __init__(self, keep: str) -> None:
        self._keep = keep

    def accept_chain(self, chain: Chain) -> bool:
        return chain.id == self._keep


def _write_subset(structure: Structure, selector: Select, path: Path) -> None:
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(path), select=selector)


def _sasa(pdb_path: Path) -> float:
    structure = freesasa.Structure(str(pdb_path))
    result = freesasa.calc(structure)
    return float(result.totalArea())


def compute_bsa(complex_pdb: Path) -> float:
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(complex_pdb.stem, str(complex_pdb))
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        pmhc_pdb = tmp_dir / "pmhc.pdb"
        binder_pdb = tmp_dir / "binder.pdb"
        _write_subset(structure, _DropChain("D"), pmhc_pdb)
        _write_subset(structure, _KeepChain("D"), binder_pdb)
        sasa_complex = _sasa(complex_pdb)
        sasa_pmhc = _sasa(pmhc_pdb)
        sasa_binder = _sasa(binder_pdb)
    return sasa_pmhc + sasa_binder - sasa_complex


def find_artifact(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if not matches:
        matches = sorted(Path(p) for p in glob(str(directory / pattern)))
    if not matches:
        raise FileNotFoundError(f"no file matching {pattern!r} under {directory}")
    return matches[0]


def chain_ids_per_residue(pdb_path: Path) -> np.ndarray:
    """Return an array of chain ids in PDB residue iteration order.

    Matches the row/column ordering of ColabFold's PAE matrix and the
    pLDDT array — both follow Bio.PDB iteration order on the rank_001
    PDB. HETATM/water residues are dropped.
    """
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = structure[0]
    ids: list[str] = []
    for chain in model.get_chains():
        for residue in chain.get_residues():
            if residue.id[0].strip() == "":
                ids.append(chain.id)
    return np.asarray(ids, dtype="<U1")


def load_plddt_from_pdb(pdb_path: Path) -> np.ndarray:
    """Per-residue pLDDT from the CA B-factor column of an AF2 PDB.

    AF2 / ColabFold encode pLDDT in the B-factor field of every atom of
    a residue; we read the CA atom. Residues lacking CA fall back to the
    first heavy atom's B-factor, which keeps the array length aligned
    with the PAE matrix even for pathological inputs.
    """
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = structure[0]
    values: list[float] = []
    for chain in model.get_chains():
        for residue in chain.get_residues():
            if residue.id[0].strip() != "":
                continue
            if "CA" in residue:
                values.append(float(residue["CA"].get_bfactor()))
                continue
            first = next(iter(residue.get_atoms()), None)
            values.append(float(first.get_bfactor()) if first is not None else 0.0)
    return np.asarray(values, dtype=np.float64)


def compute_heavy_atom_distances(pdb_path: Path) -> np.ndarray:
    """Pairwise residue-residue min heavy-atom distance, shape (N, N).

    Residue order matches ``chain_ids_per_residue``. Hydrogens are
    excluded via ``atom.element != "H"``. Diagonal is 0; residues with
    no heavy atoms keep +inf rows/cols.
    """
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = structure[0]
    residue_coords: list[np.ndarray] = []
    for chain in model.get_chains():
        for residue in chain.get_residues():
            if residue.id[0].strip() != "":
                continue
            coords = [a.get_coord() for a in residue.get_atoms() if a.element != "H"]
            if coords:
                residue_coords.append(np.asarray(coords, dtype=np.float64))
            else:
                residue_coords.append(np.zeros((0, 3), dtype=np.float64))
    n = len(residue_coords)
    distances = np.full((n, n), np.inf, dtype=np.float64)
    for i in range(n):
        ai = residue_coords[i]
        if ai.shape[0] == 0:
            continue
        for j in range(i, n):
            aj = residue_coords[j]
            if aj.shape[0] == 0:
                continue
            d = float(np.linalg.norm(ai[:, None, :] - aj[None, :, :], axis=-1).min())
            distances[i, j] = d
            distances[j, i] = d
    np.fill_diagonal(distances, 0.0)
    return distances


def _interface_mask(
    chain_ids: np.ndarray,
    distances: np.ndarray,
    binder: str,
    target: Sequence[str],
    thr: float,
) -> np.ndarray:
    binder_mask = chain_ids == binder
    target_mask = np.isin(chain_ids, list(target))
    cross = binder_mask[:, None] & target_mask[None, :]
    return cross & (distances <= thr)


def compute_ipae_interface(
    pae: list[list[float]] | np.ndarray,
    chain_ids: np.ndarray,
    distances: np.ndarray,
    binder: str = "D",
    target: Sequence[str] = ("A", "B", "C"),
    thr: float = 8.0,
) -> float:
    """Canonical iPAE: mean min(pae_ij, pae_ji) over D <-> A/B/C residue
    pairs whose heavy-atom min distance is <= ``thr`` (8 A default).

    Returns +inf when no interface exists — signals a failed design.
    """
    pae_arr = np.asarray(pae, dtype=np.float64)
    mask = _interface_mask(chain_ids, distances, binder, target, thr)
    if not mask.any():
        return float("inf")
    combined = np.minimum(pae_arr, pae_arr.T)
    return float(np.mean(combined[mask]))


def compute_iplddt_interface(
    plddt: list[float] | np.ndarray,
    chain_ids: np.ndarray,
    distances: np.ndarray,
    binder: str = "D",
    target: Sequence[str] = ("A", "B", "C"),
    thr: float = 8.0,
) -> float:
    """Mean pLDDT over residues that touch the binder<->target interface
    (rows OR columns of the interface mask). Returns 0.0 when no
    interface exists — pairs with iPAE=+inf to mark a failed design.
    """
    plddt_arr = np.asarray(plddt, dtype=np.float64)
    mask = _interface_mask(chain_ids, distances, binder, target, thr)
    if not mask.any():
        return 0.0
    touched = mask.any(axis=1) | mask.any(axis=0)
    return float(np.mean(plddt_arr[touched]))


def metrics_for_design(pred_dir: Path) -> dict[str, Any]:
    """Canonical metrics for one Stage 2 designs AF2 prediction.

    Reads the rank_001 PDB (chain ids, heavy-atom geometry, CA B-factor
    pLDDT) and the scores JSON (PAE matrix); routes through the
    interface-mask iPAE/ipLDDT. BSA computation is shared with the
    controls path.
    """
    scores_json = find_artifact(pred_dir, DEFAULT_SCORES_GLOB)
    rank_pdb = find_artifact(pred_dir, DEFAULT_PDB_GLOB)
    pae = load_pae(scores_json)
    chain_ids = chain_ids_per_residue(rank_pdb)
    distances = compute_heavy_atom_distances(rank_pdb)
    plddt = load_plddt_from_pdb(rank_pdb)
    ipae = compute_ipae_interface(pae, chain_ids, distances)
    iplddt = compute_iplddt_interface(plddt, chain_ids, distances)
    bsa = compute_bsa(rank_pdb)
    return {
        "ipae": ipae,
        "iplddt": iplddt,
        "bsa": bsa,
        "rank_001_pdb": str(rank_pdb),
        "scores_json": str(scores_json),
    }


def metrics_for_control(colabfold_dir: Path, backbone_pdb: Path | None) -> dict[str, Any]:
    scores_json = find_artifact(colabfold_dir, DEFAULT_SCORES_GLOB)
    rank_pdb = find_artifact(colabfold_dir, DEFAULT_PDB_GLOB)
    # Chain lengths come from the predicted complex PDB itself — the cleaned
    # crystal structure (when provided) lacks chain D, so we always slice from
    # the ColabFold rank_001 PDB, which has the full A+B+C+D layout.
    lengths = chain_lengths(rank_pdb)
    slices = chain_slices(lengths)
    pae = load_pae(scores_json)
    plddt = load_plddt(scores_json)
    ipae = compute_ipae(pae, slices)
    iplddt = compute_iplddt(plddt, slices)
    bsa = compute_bsa(rank_pdb)
    return {
        "ipae": ipae,
        "iplddt": iplddt,
        "bsa": bsa,
        "rank_001_pdb": str(rank_pdb),
        "scores_json": str(scores_json),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--colabfold-dir", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--mock", action="store_true")
    # Back-compat with original stub signature.
    parser.add_argument("--config")
    args = parser.parse_args(argv)

    colabfold_dir = args.colabfold_dir or (Path(args.config).parent if args.config else None)
    out_path = args.out
    if colabfold_dir is None or out_path is None:
        parser.error("must provide --colabfold-dir and --out (or legacy --config/--out)")

    if args.mock:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(
            json.dumps({"mock": True, "controls": {}, "designs": []}, indent=2) + "\n"
        )
        print(f"[mock] compute_metrics: wrote stub {out_path}")
        return 0

    result = metrics_for_control(colabfold_dir, None)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
