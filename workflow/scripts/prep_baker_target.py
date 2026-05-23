# mypy: disable-error-code="no-untyped-call,attr-defined,no-any-return,misc"
"""Cycle 03 sub-run A -- BAKER-format truncated target.

BAKER lab's published pMHC-I binder scaffolds (BAKER_LAB_2025) carry the
target as a single chain B = HLA-A*02:01 alpha1+alpha2 (residues 1-180)
fused directly to the 9-mer peptide -- no beta2-microglobulin, no alpha3
domain. Our canonical ``3hpj_clean.pdb`` instead follows the crystal-frame
convention (chain A = HC ~275 res, B = beta2m, C = peptide). Aligning a
BAKER scaffold's HLA-bearing chain onto our beta2m chain silently superposes
two unrelated folds, so sub-run A must align against a target whose chain
layout matches BAKER's groove-only fragment.

This script produces that fragment from ``3hpj_clean.pdb``:

* keep chain A residues 1-180 (the alpha1+alpha2 peptide-binding groove),
  drop 181-275 (alpha3),
* rename that chain A -> B (BAKER convention),
* drop the original chain B (beta2m) entirely,
* keep chain C (peptide, residues 1-9) verbatim,
* ATOM records only (HETATM/water/H/alt-loc-B stripped), original numbering
  preserved (no renumber).

Output: ``data/targets/3hpj_baker_truncated.pdb`` with exactly 180 CA on
chain B and 9 CA on chain C. Sub-run B keeps the full target unchanged.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
from pathlib import Path

from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from workflow.scripts.prep_target import _Keep, _rename_chains_safely, _residues

logger = logging.getLogger("prep_baker_target")
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

DEFAULT_INPUT = Path("data/targets/3hpj_clean.pdb")
DEFAULT_OUTPUT = Path("data/targets/3hpj_baker_truncated.pdb")
MOCK_FIXTURE = Path("tests/fixtures/targets/3hpj_baker_truncated_mock.pdb")

HLA_GROOVE_LAST_RESNUM = 180  # alpha1+alpha2; alpha3 starts at 181
PEPTIDE_LENGTH = 9


def _truncate_chain_a(model: Model, last_resnum: int) -> None:
    """Drop chain A residues numbered above ``last_resnum`` (the alpha3 domain)
    and any HETATM, in place. Standard residues 1..last_resnum are kept."""
    chain = model["A"]
    for residue in list(chain.get_residues()):
        hetflag, resnum, _ = residue.id
        if hetflag.strip() != "" or resnum > last_resnum:
            chain.detach_child(residue.id)


def truncate_target(
    clean_pdb: Path,
    out_pdb: Path,
    last_resnum: int = HLA_GROOVE_LAST_RESNUM,
    peptide_length: int = PEPTIDE_LENGTH,
) -> dict[str, int]:
    """Truncate the canonical clean target into the BAKER groove-only layout.

    Returns the per-chain CA count of the written PDB ({"B": 180, "C": 9}).
    """
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure("3hpj_baker_truncated", str(clean_pdb))
    model: Model = structure[0]

    present = {ch.id for ch in model.get_chains()}
    missing = [c for c in ("A", "B", "C") if c not in present]
    if missing:
        raise ValueError(f"{clean_pdb}: chains {missing} absent; expected canonical A/B/C")

    model.detach_child("B")  # drop beta2m
    _truncate_chain_a(model, last_resnum)

    hla_len = len(_residues(model["A"]))
    pep_len = len(_residues(model["C"]))
    if hla_len != last_resnum:
        raise ValueError(
            f"{clean_pdb}: HLA groove chain A has {hla_len} residues 1..{last_resnum}, "
            f"expected {last_resnum} (renumbered or missing residues?)"
        )
    if pep_len != peptide_length:
        raise ValueError(
            f"{clean_pdb}: peptide chain C length {pep_len} != expected {peptide_length}"
        )

    _rename_chains_safely(model, {"A": "B"})  # HLA groove A -> B (BAKER convention)

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_pdb), select=_Keep({"B", "C"}))

    counts = _ca_counts(out_pdb)
    expected = {"B": last_resnum, "C": peptide_length}
    if counts != expected:
        raise ValueError(f"{out_pdb}: CA chain counts {counts} != expected {expected}")
    logger.info("wrote %s (CA counts %s)", out_pdb, counts)
    return counts


def _ca_counts(pdb_path: Path) -> dict[str, int]:
    parser = PDBParser(QUIET=True)
    model = parser.get_structure(pdb_path.stem, str(pdb_path))[0]
    counts: dict[str, int] = {}
    for chain in model.get_chains():
        n = sum(1 for r in chain.get_residues() if r.id[0].strip() == "" and "CA" in r)
        counts[chain.id] = n
    return counts


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in", dest="clean_pdb", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out", dest="out_pdb", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--hla-last-resnum", type=int, default=HLA_GROOVE_LAST_RESNUM)
    parser.add_argument("--peptide-length", type=int, default=PEPTIDE_LENGTH)
    parser.add_argument("--mock", action="store_true")
    args = parser.parse_args(argv)

    if args.mock:
        args.out_pdb.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(MOCK_FIXTURE, args.out_pdb)
        logger.info("mock: copied %s -> %s", MOCK_FIXTURE, args.out_pdb)
        return 0

    truncate_target(args.clean_pdb, args.out_pdb, args.hla_last_resnum, args.peptide_length)
    return 0


if __name__ == "__main__":
    sys.exit(main())
