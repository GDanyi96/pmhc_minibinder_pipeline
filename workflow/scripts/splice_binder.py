# mypy: disable-error-code="no-untyped-call,attr-defined,no-any-return,misc"
"""Stage 2 designs — splice a Stage 1 binder onto the cleaned pMHC.

Stage 1 (RFdiffusion) emits each binder as chain "A" inside its own
single-chain PDB. Stage 2 ProteinMPNN expects a 4-chain complex where the
fixed scaffold (A=HC, B=beta2m, C=peptide) is the cleaned pMHC and the
designed chain is D=binder. This module reads both PDBs, renames the
Stage 1 chain "A" to "D", composes the four chains into one structure
with per-chain residue renumbering starting at 1, and writes the spliced
PDB to disk.

Reference: docs/known_traps.md trap #17 — skipping the rename produces a
3-chain AF2 input where the binder is silently treated as HC, yielding
garbage iPAE.
"""

from __future__ import annotations

from pathlib import Path

from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure

_STAGE1_BINDER_CHAIN = "A"
_TARGET_BINDER_CHAIN = "D"
_PMHC_CHAINS = ("A", "B", "C")


def _standard_residues(chain: Chain) -> list:  # type: ignore[type-arg]
    """Residues with empty HETATM flag (drops waters, ions, etc.)."""
    return [r for r in chain.get_residues() if r.id[0].strip() == ""]


def _load_first_model(pdb_path: Path) -> Model:
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    return structure[0]


def _renumber_chain(chain: Chain) -> None:
    """Renumber standard residues 1..N in place; drop HETATMs."""
    residues = _standard_residues(chain)
    for r in list(chain.get_residues()):
        chain.detach_child(r.id)
    for new_resnum, residue in enumerate(residues, start=1):
        residue.id = (" ", new_resnum, " ")
        chain.add(residue)


def splice_binder_onto_pmhc(
    stage1_pdb: Path,
    cleaned_pmhc_pdb: Path,
    out_pdb: Path,
) -> None:
    """Compose a 4-chain (A=HC, B=beta2m, C=peptide, D=binder) PDB.

    Asserts the stage1 PDB has exactly one chain (the binder) and the
    cleaned pMHC has exactly A/B/C. Output chains are renumbered 1..N
    per chain — AF2 multimer input contract.
    """
    s1_model = _load_first_model(stage1_pdb)
    s1_chains = list(s1_model.get_chains())
    if len(s1_chains) != 1:
        raise ValueError(
            f"{stage1_pdb}: expected exactly 1 chain (the binder), found "
            f"{[c.id for c in s1_chains]}"
        )
    if s1_chains[0].id != _STAGE1_BINDER_CHAIN:
        raise ValueError(
            f"{stage1_pdb}: expected binder on chain {_STAGE1_BINDER_CHAIN!r}, "
            f"found chain {s1_chains[0].id!r}"
        )

    pmhc_model = _load_first_model(cleaned_pmhc_pdb)
    pmhc_chain_ids = {c.id for c in pmhc_model.get_chains()}
    if _TARGET_BINDER_CHAIN in pmhc_chain_ids:
        raise ValueError(
            f"{cleaned_pmhc_pdb}: chain {_TARGET_BINDER_CHAIN!r} already present; "
            "cleaned pMHC must contain only A/B/C (HC, beta2m, peptide)"
        )
    missing = [c for c in _PMHC_CHAINS if c not in pmhc_chain_ids]
    if missing:
        raise ValueError(
            f"{cleaned_pmhc_pdb}: missing required chains {missing}; "
            f"found {sorted(pmhc_chain_ids)}"
        )

    binder_chain = s1_chains[0]
    binder_chain.id = _TARGET_BINDER_CHAIN

    out_structure = Structure(out_pdb.stem)
    out_model = Model(0)
    out_structure.add(out_model)
    for cid in _PMHC_CHAINS:
        chain = pmhc_model[cid].copy()
        _renumber_chain(chain)
        out_model.add(chain)
    binder_copy = binder_chain.copy()
    _renumber_chain(binder_copy)
    out_model.add(binder_copy)

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(out_structure)
    io.save(str(out_pdb))
