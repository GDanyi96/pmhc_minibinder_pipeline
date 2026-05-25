# mypy: disable-error-code="no-untyped-call,attr-defined,no-any-return,misc"
"""Cycle 03 -- BioPython peptide-contact filter.

Cycle 02's placement deficit (65% of designs made zero hotspot contact)
motivates an explicit, Rosetta-free contact gate between ProteinMPNN and
AF2 in the cycle-03 funnel: a design is kept only if its binder (chain D)
makes at least ``min_contacts`` C-beta contacts within ``distance`` A of any
peptide (chain C) atom. This cheaply discards mis-placed binders before the
expensive AF2 fan-in, rather than discovering the failure post-fold.

Glycine has no C-beta; we fall back to its C-alpha so every binder residue
contributes one representative atom. Thresholds live in
configs/thresholds.yaml (no magic numbers) and are passed in by the caller;
the module-level constants are documentation defaults only.

This is a standalone prep component: it is unit-tested and spec'd here, and
wired into the real cycle-03 Stage 2 funnel (MPNN -> contact_filter -> AF2);
it is intentionally NOT inserted into the mock Snakemake DAG so the existing
Stage 2 4/40 halt-gate calibration is preserved.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure

logger = logging.getLogger("contact_filter")
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

DEFAULT_DISTANCE = 5.0
DEFAULT_MIN_CONTACTS = 3
_BINDER_CHAIN = "D"
_PEPTIDE_CHAIN = "C"


def _representative_coords(structure: Structure, chain_id: str) -> np.ndarray:
    """C-beta coords for each standard residue of ``chain_id`` (C-alpha for Gly)."""
    model = structure[0]
    if chain_id not in model:
        raise ValueError(f"chain {chain_id!r} absent from {structure.id}")
    coords: list[list[float]] = []
    for residue in model[chain_id].get_residues():
        if residue.id[0].strip() != "":
            continue
        atom = None
        if "CB" in residue:
            atom = residue["CB"]
        elif "CA" in residue:
            atom = residue["CA"]
        if atom is not None:
            coords.append([float(c) for c in atom.get_coord()])
    return np.asarray(coords, dtype=np.float64)


def _all_atom_coords(structure: Structure, chain_id: str) -> np.ndarray:
    model = structure[0]
    if chain_id not in model:
        raise ValueError(f"chain {chain_id!r} absent from {structure.id}")
    coords = [
        [float(c) for c in atom.get_coord()]
        for residue in model[chain_id].get_residues()
        if residue.id[0].strip() == ""
        for atom in residue.get_atoms()
    ]
    return np.asarray(coords, dtype=np.float64)


def count_peptide_contacts(
    pdb_path: Path,
    distance: float = DEFAULT_DISTANCE,
    binder_chain: str = _BINDER_CHAIN,
    peptide_chain: str = _PEPTIDE_CHAIN,
) -> int:
    """Number of binder C-beta atoms within ``distance`` A of any peptide atom."""
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    binder_cb = _representative_coords(structure, binder_chain)
    peptide_atoms = _all_atom_coords(structure, peptide_chain)
    if binder_cb.size == 0 or peptide_atoms.size == 0:
        return 0
    dists = np.linalg.norm(binder_cb[:, None, :] - peptide_atoms[None, :, :], axis=-1)
    return int((dists.min(axis=1) <= distance).sum())


def passes(
    pdb_path: Path,
    distance: float = DEFAULT_DISTANCE,
    min_contacts: int = DEFAULT_MIN_CONTACTS,
    binder_chain: str = _BINDER_CHAIN,
    peptide_chain: str = _PEPTIDE_CHAIN,
) -> bool:
    return count_peptide_contacts(pdb_path, distance, binder_chain, peptide_chain) >= min_contacts


def filter_designs(
    pdbs: list[Path],
    distance: float = DEFAULT_DISTANCE,
    min_contacts: int = DEFAULT_MIN_CONTACTS,
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for pdb in pdbs:
        n = count_peptide_contacts(pdb, distance)
        records.append(
            {
                "pdb": str(pdb),
                "peptide_contacts": n,
                "pass": n >= min_contacts,
            }
        )
    return records


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-glob", default=None)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--distance", type=float, default=DEFAULT_DISTANCE)
    parser.add_argument("--min-contacts", type=int, default=DEFAULT_MIN_CONTACTS)
    parser.add_argument("--mock", action="store_true")
    args = parser.parse_args(argv)

    glob = args.input_glob
    if glob is None:
        # The filter operates on 4-chain spliced designs (A/B/C/D); the mock
        # default is the stitched seed complex fixture, real default is the
        # Stage 2 spliced output.
        glob = (
            "tests/fixtures/design_2079_mock_seed.pdb"
            if args.mock
            else "results/**/spliced/design_*.pdb"
        )
    base = Path(glob)
    pdbs = sorted(base.parent.glob(base.name))
    records = filter_designs(pdbs, args.distance, args.min_contacts)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")
    n_pass = sum(1 for r in records if r["pass"])
    logger.info("contact filter: %d/%d designs pass", n_pass, len(records))
    return 0


if __name__ == "__main__":
    sys.exit(main())
