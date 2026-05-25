# mypy: disable-error-code="no-untyped-call,attr-defined,no-any-return,misc"
"""Stage 2 designs — splice a Stage 1 binder onto the cleaned pMHC.

Stage 1 (RFdiffusion) emits a 4-chain complex (A=HC, B=beta2m, C=peptide,
D=binder) under the cycle 02 contig
``[A1-275/0 B0-99/0 C1-9/0 70-110]``. The motif chains A/B/C are passed
through unchanged by RFdiffusion and should match the canonical cleaned
pMHC by construction (motif RMSD ~0.12 Å verified empirically on cycle
02). Stage 2 ProteinMPNN expects a 4-chain complex where the fixed
scaffold (A=HC, B=beta2m, C=peptide) is the *canonical* cleaned pMHC and
the designed chain is D=binder. This module extracts chain D from the
Stage 1 PDB, combines it with A/B/C from the cleaned pMHC, renumbers
each chain 1..N, and writes the composed PDB to disk.

References: docs/known_traps.md trap #17 (without chain D the AF2
input is silently treated as 3-chain HC+B2M+peptide, yielding garbage
iPAE); trap #28 (RFdiffusion writes binder as the next-free letter, not
chain A — false belief previously embedded in this module); trap #29
(the four-chain rewrite below).
"""

from __future__ import annotations

from pathlib import Path

from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure

_STAGE1_BINDER_CHAIN = "D"
_TARGET_BINDER_CHAIN = "D"
_PMHC_CHAINS = ("A", "B", "C")
_ALLOWED_STAGE1_CHAINS = frozenset({"A", "B", "C", "D"})

# Cycle 03 sub-run A (truncated, 3-chain). RFdiffusion partial diffusion
# preserves the aligned scaffold's chain IDs: A=binder, B=HLA[1:180], C=peptide.
# The truncated cleaned target carries HLA on B and peptide on C (no chain A,
# no beta2m). See PROJECT_STATE 12.j and Trap #41.
_SUBRUN_A_DESIGN_BINDER_CHAIN = "A"
_ALLOWED_SUBRUN_A_CHAINS = frozenset({"A", "B", "C"})
_TRUNCATED_TARGET_HLA_CHAIN = "B"
_TRUNCATED_TARGET_PEPTIDE_CHAIN = "C"


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

    Asserts the stage1 PDB contains a chain D (the binder) and no
    chains outside {A,B,C,D}; the chains A/B/C present in the stage1
    PDB are ignored in favour of the canonical cleaned pMHC. The
    cleaned pMHC must have exactly A/B/C (and no D). Output chains are
    renumbered 1..N per chain — AF2 multimer input contract.
    """
    s1_model = _load_first_model(stage1_pdb)
    s1_chain_ids = {c.id for c in s1_model.get_chains()}
    if _STAGE1_BINDER_CHAIN not in s1_chain_ids:
        raise ValueError(
            f"{stage1_pdb}: expected binder on chain {_STAGE1_BINDER_CHAIN!r}, "
            f"found chains {sorted(s1_chain_ids)}"
        )
    unexpected = s1_chain_ids - _ALLOWED_STAGE1_CHAINS
    if unexpected:
        raise ValueError(
            f"{stage1_pdb}: unexpected chain(s) {sorted(unexpected)}; "
            f"allowed chains are {sorted(_ALLOWED_STAGE1_CHAINS)}"
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

    binder_chain = s1_model[_STAGE1_BINDER_CHAIN]

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


def splice_binder_subrun_a(
    design_pdb: Path,
    cleaned_truncated_pdb: Path,
    out_pdb: Path,
) -> None:
    """Compose a 3-chain truncated complex A=HLA, B=peptide, C=binder.

    Cycle 03 sub-run A counterpart to ``splice_binder_onto_pmhc``. The sub-run A
    RFdiffusion output is 3-chain (A=binder, B=HLA[1:180], C=peptide); the
    truncated cleaned target carries HLA on chain B and peptide on chain C (no
    beta2m, no chain A). We take the binder from the design's chain A and the
    HLA/peptide from the *canonical* truncated target (coordinate frame is
    preserved through alignment + partial diffusion, motif RMSD ~0.11 A), then
    emit them in AF2 FASTA order A=HLA, B=peptide, C=binder so the ColabFold
    output matches ``LAYOUT_CHAINS["truncated"]`` (binder = C). Output chains are
    renumbered 1..N per chain.
    """
    design_model = _load_first_model(design_pdb)
    design_chain_ids = {c.id for c in design_model.get_chains()}
    if _SUBRUN_A_DESIGN_BINDER_CHAIN not in design_chain_ids:
        raise ValueError(
            f"{design_pdb}: expected sub-run A binder on chain "
            f"{_SUBRUN_A_DESIGN_BINDER_CHAIN!r}, found chains {sorted(design_chain_ids)}"
        )
    unexpected = design_chain_ids - _ALLOWED_SUBRUN_A_CHAINS
    if unexpected:
        raise ValueError(
            f"{design_pdb}: unexpected chain(s) {sorted(unexpected)}; "
            f"allowed sub-run A chains are {sorted(_ALLOWED_SUBRUN_A_CHAINS)}"
        )

    pmhc_model = _load_first_model(cleaned_truncated_pdb)
    pmhc_chain_ids = {c.id for c in pmhc_model.get_chains()}
    missing = [
        c
        for c in (_TRUNCATED_TARGET_HLA_CHAIN, _TRUNCATED_TARGET_PEPTIDE_CHAIN)
        if c not in pmhc_chain_ids
    ]
    if missing:
        raise ValueError(
            f"{cleaned_truncated_pdb}: missing required truncated target chain(s) "
            f"{missing}; found {sorted(pmhc_chain_ids)} "
            "(expected B=HLA[1:180], C=peptide)"
        )

    binder_chain = design_model[_SUBRUN_A_DESIGN_BINDER_CHAIN]

    out_structure = Structure(out_pdb.stem)
    out_model = Model(0)
    out_structure.add(out_model)
    # source chain id -> output chain id
    composition = [
        (pmhc_model[_TRUNCATED_TARGET_HLA_CHAIN], "A"),
        (pmhc_model[_TRUNCATED_TARGET_PEPTIDE_CHAIN], "B"),
        (binder_chain, "C"),
    ]
    for src_chain, out_cid in composition:
        chain = src_chain.copy()
        chain.id = out_cid
        _renumber_chain(chain)
        out_model.add(chain)

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(out_structure)
    io.save(str(out_pdb))
