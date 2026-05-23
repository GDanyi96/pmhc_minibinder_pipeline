# mypy: disable-error-code="no-untyped-call,attr-defined,no-any-return,misc"
"""Cycle 03 sub-run A -- align BAKER scaffolds onto our target frame.

BAKER lab's published pMHC-I binder scaffolds (BAKER_LAB_2025) live in their
own crystallographic frame. Before partial diffusion (partial_diffuse.py) we
rigid-body align each scaffold's target chain onto our cleaned reference
``3hpj_clean.pdb`` so the diffused binder lands against the WT1/A*02:01
interface in our coordinate system.

This is a BioPython ``Superimposer`` equivalent of BAKER's ``align_chainB.py``
(``/workspace/pMHCI_binder_design/``); we avoid the upstream-script dependency
for the same reason we avoid Rosetta -- it keeps the stage runnable anywhere,
mock or pod, with no external clone on PATH. The superposition matches CA
atoms of the scaffold target chain to the reference MHC chain (first
``min(len)`` CA atoms, in order), then applies the recovered rotation +
translation to every atom of the scaffold.

Mock and real mode share the same code; only the input glob differs
(tests/fixtures/baker_library_mock/ vs data/scaffolds/baker_library/).
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from Bio.PDB import PDBIO, PDBParser, Superimposer
from Bio.PDB.Structure import Structure

from workflow.scripts import align_baker_scaffolds

logger = logging.getLogger("align_scaffolds")
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

DEFAULT_SCAFFOLD_TARGET_CHAIN = "B"
DEFAULT_REFERENCE_CHAIN = "A"


def _ca_atoms(structure: Structure, chain_id: str) -> list:  # type: ignore[type-arg]
    model = structure[0]
    if chain_id not in model:
        raise ValueError(f"chain {chain_id!r} absent from {structure.id}")
    return [r["CA"] for r in model[chain_id].get_residues() if r.id[0].strip() == "" and "CA" in r]


def align_one(
    scaffold_pdb: Path,
    reference_pdb: Path,
    out_pdb: Path,
    scaffold_target_chain: str = DEFAULT_SCAFFOLD_TARGET_CHAIN,
    reference_chain: str = DEFAULT_REFERENCE_CHAIN,
) -> float:
    """Align one scaffold onto the reference; write it and return the RMSD."""
    parser = PDBParser(QUIET=True)
    scaffold = parser.get_structure(scaffold_pdb.stem, str(scaffold_pdb))
    reference = parser.get_structure(reference_pdb.stem, str(reference_pdb))

    moving = _ca_atoms(scaffold, scaffold_target_chain)
    fixed = _ca_atoms(reference, reference_chain)
    k = min(len(moving), len(fixed))
    if k < 3:
        raise ValueError(
            f"{scaffold_pdb}: need >=3 matched CA atoms to superpose "
            f"(scaffold chain {scaffold_target_chain}={len(moving)}, "
            f"reference chain {reference_chain}={len(fixed)})"
        )

    sup = Superimposer()
    sup.set_atoms(fixed[:k], moving[:k])
    sup.apply(scaffold.get_atoms())

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(scaffold)
    io.save(str(out_pdb))
    return 0.0 if sup.rms is None else float(sup.rms)


def align_scaffolds(
    scaffold_glob: str,
    reference_pdb: Path,
    out_dir: Path,
    scaffold_target_chain: str = DEFAULT_SCAFFOLD_TARGET_CHAIN,
    reference_chain: str = DEFAULT_REFERENCE_CHAIN,
    baker_layout: bool = False,
) -> list[Path]:
    """Align every scaffold matching ``scaffold_glob`` onto ``reference_pdb``.

    Returns the aligned scaffold paths (sorted, deterministic order).

    ``baker_layout=True`` delegates to align_baker_scaffolds.py: BAKER's
    published library carries the target as a fused chain B (HLA[1:180] +
    peptide) with the binder on chain A, so it must be aligned against the
    truncated reference (chain B = HLA, chain C = peptide) and rewritten into
    our A=binder / B=HLA / C=peptide layout. Used by cycle-03 sub-run A.
    """
    if baker_layout:
        return align_baker_scaffolds.align_baker_scaffolds(scaffold_glob, reference_pdb, out_dir)
    base = Path(scaffold_glob)
    matches = sorted(base.parent.glob(base.name))
    if not matches:
        raise FileNotFoundError(f"no scaffolds matched {scaffold_glob}")
    out_dir.mkdir(parents=True, exist_ok=True)
    aligned: list[Path] = []
    for scaffold in matches:
        out_pdb = out_dir / scaffold.name
        rms = align_one(scaffold, reference_pdb, out_pdb, scaffold_target_chain, reference_chain)
        logger.info("aligned %s -> %s (CA RMSD %.3f A)", scaffold.name, out_pdb, rms)
        aligned.append(out_pdb)
    return aligned


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scaffold-glob", default=None)
    parser.add_argument("--reference-pdb", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--scaffold-target-chain", default=DEFAULT_SCAFFOLD_TARGET_CHAIN)
    parser.add_argument("--reference-chain", default=DEFAULT_REFERENCE_CHAIN)
    parser.add_argument("--mock", action="store_true")
    args = parser.parse_args(argv)

    if args.mock:
        scaffold_glob = args.scaffold_glob or "tests/fixtures/baker_library_mock/scaf*.pdb"
        reference_pdb = args.reference_pdb or Path("tests/fixtures/stage1/mock_clean.pdb")
    else:
        scaffold_glob = args.scaffold_glob or "data/scaffolds/baker_library/scaf*.pdb"
        reference_pdb = args.reference_pdb or Path("data/targets/3hpj_clean.pdb")

    aligned = align_scaffolds(
        scaffold_glob,
        reference_pdb,
        args.out_dir,
        args.scaffold_target_chain,
        args.reference_chain,
    )
    logger.info("aligned %d scaffolds into %s", len(aligned), args.out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
