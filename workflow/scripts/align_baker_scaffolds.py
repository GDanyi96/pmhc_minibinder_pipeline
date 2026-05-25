# mypy: disable-error-code="no-untyped-call,attr-defined,no-any-return,misc"
"""Cycle 03 sub-run A -- align BAKER scaffolds (fused chain B) onto our target.

BAKER lab's published pMHC-I binder scaffolds (BAKER_LAB_2025) carry the
target as a single chain B = HLA-A*02:01 alpha1+alpha2 (residues 1-180)
fused directly to the 9-mer peptide (residues 181-189); chain A is the
binder. Our reference is ``3hpj_baker_truncated.pdb`` (chain B = HLA[1:180],
chain C = peptide[1:9]; see prep_baker_target.py).

For each scaffold we superpose the scaffold's chain-B HLA residues (1-180)
onto the reference chain-B HLA residues (1-180) -- the first 180 CA, in
order, which are sequence-identical -- then apply the recovered rigid
transform to the whole scaffold so the binder lands in our target frame.
The output is rewritten into our standard layout: chain A = transformed
binder, chain B = reference HLA (180 res), chain C = reference peptide
(9 res). BAKER's fused chain B and any extra context are discarded in favour
of our canonical reference chains.

This is the BAKER-fused-layout counterpart of align_scaffolds.py and a
BioPython ``Superimposer`` equivalent of BAKER's upstream ``align_chainB.py``
(``/workspace/pMHCI_binder_design/``); avoiding the external clone keeps the
stage runnable anywhere, mock or pod, as with align_scaffolds.py.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from Bio.PDB import PDBIO, PDBParser, Superimposer
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure

logger = logging.getLogger("align_baker_scaffolds")
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

DEFAULT_SCAFFOLD_BINDER_CHAIN = "A"
DEFAULT_SCAFFOLD_FUSED_CHAIN = "B"
DEFAULT_REFERENCE_HLA_CHAIN = "B"
DEFAULT_REFERENCE_PEPTIDE_CHAIN = "C"


def _chain_ca(model: Model, chain_id: str) -> list:  # type: ignore[type-arg]
    if chain_id not in model:
        raise ValueError(f"chain {chain_id!r} absent from {model.get_parent().id}")
    return [r["CA"] for r in model[chain_id].get_residues() if r.id[0].strip() == "" and "CA" in r]


def align_one_baker(
    scaffold_pdb: Path,
    reference_pdb: Path,
    out_pdb: Path,
    scaffold_binder_chain: str = DEFAULT_SCAFFOLD_BINDER_CHAIN,
    scaffold_fused_chain: str = DEFAULT_SCAFFOLD_FUSED_CHAIN,
    reference_hla_chain: str = DEFAULT_REFERENCE_HLA_CHAIN,
    reference_peptide_chain: str = DEFAULT_REFERENCE_PEPTIDE_CHAIN,
) -> float:
    """Align one BAKER scaffold onto the truncated reference; write the
    A=binder / B=HLA / C=peptide complex and return the superposition RMSD."""
    parser = PDBParser(QUIET=True)
    scaffold = parser.get_structure(scaffold_pdb.stem, str(scaffold_pdb))
    reference = parser.get_structure(reference_pdb.stem, str(reference_pdb))
    s_model: Model = scaffold[0]
    r_model: Model = reference[0]

    # Reference HLA is 180 res; the scaffold's fused chain B is HLA(180) +
    # peptide(9) = 189. min() selects the first 180 (the HLA substring), which
    # is sequence-identical to the reference HLA, so position-order CA matching
    # is the correct substring superposition.
    moving = _chain_ca(s_model, scaffold_fused_chain)
    fixed = _chain_ca(r_model, reference_hla_chain)
    k = min(len(moving), len(fixed))
    if k < 3:
        raise ValueError(
            f"{scaffold_pdb}: need >=3 matched HLA CA atoms to superpose "
            f"(scaffold chain {scaffold_fused_chain}={len(moving)}, "
            f"reference chain {reference_hla_chain}={len(fixed)})"
        )

    sup = Superimposer()
    sup.set_atoms(fixed[:k], moving[:k])
    sup.apply(s_model.get_atoms())

    if scaffold_binder_chain not in s_model:
        raise ValueError(f"{scaffold_pdb}: binder chain {scaffold_binder_chain!r} absent")

    out_structure = Structure(out_pdb.stem)
    out_model = Model(0)
    out_structure.add(out_model)
    binder = s_model[scaffold_binder_chain].copy()
    binder.id = "A"
    hla = r_model[reference_hla_chain].copy()
    hla.id = "B"
    peptide = r_model[reference_peptide_chain].copy()
    peptide.id = "C"
    out_model.add(binder)
    out_model.add(hla)
    out_model.add(peptide)

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(out_structure)
    io.save(str(out_pdb))
    return 0.0 if sup.rms is None else float(sup.rms)


def align_baker_scaffolds(
    scaffold_glob: str,
    reference_pdb: Path,
    out_dir: Path,
    scaffold_binder_chain: str = DEFAULT_SCAFFOLD_BINDER_CHAIN,
    scaffold_fused_chain: str = DEFAULT_SCAFFOLD_FUSED_CHAIN,
    reference_hla_chain: str = DEFAULT_REFERENCE_HLA_CHAIN,
    reference_peptide_chain: str = DEFAULT_REFERENCE_PEPTIDE_CHAIN,
) -> list[Path]:
    """Align every scaffold matching ``scaffold_glob`` onto ``reference_pdb``.

    Returns the aligned scaffold paths (sorted, deterministic order).
    """
    base = Path(scaffold_glob)
    matches = sorted(base.parent.glob(base.name))
    if not matches:
        raise FileNotFoundError(f"no scaffolds matched {scaffold_glob}")
    out_dir.mkdir(parents=True, exist_ok=True)
    aligned: list[Path] = []
    for scaffold in matches:
        out_pdb = out_dir / scaffold.name
        rms = align_one_baker(
            scaffold,
            reference_pdb,
            out_pdb,
            scaffold_binder_chain,
            scaffold_fused_chain,
            reference_hla_chain,
            reference_peptide_chain,
        )
        logger.info("aligned %s -> %s (HLA CA RMSD %.3f A)", scaffold.name, out_pdb, rms)
        aligned.append(out_pdb)
    return aligned


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scaffold-glob", default=None)
    parser.add_argument("--reference-pdb", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--mock", action="store_true")
    args = parser.parse_args(argv)

    if args.mock:
        scaffold_glob = args.scaffold_glob or "tests/fixtures/baker_truncated_mock/scaf*.pdb"
        reference_pdb = args.reference_pdb or Path(
            "tests/fixtures/targets/3hpj_baker_truncated_mock.pdb"
        )
    else:
        scaffold_glob = args.scaffold_glob or "data/scaffolds/baker_library/scaf*.pdb"
        reference_pdb = args.reference_pdb or Path("data/targets/3hpj_baker_truncated.pdb")

    aligned = align_baker_scaffolds(scaffold_glob, reference_pdb, args.out_dir)
    logger.info("aligned %d BAKER scaffolds into %s", len(aligned), args.out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
