# mypy: disable-error-code="no-untyped-call,attr-defined"
"""Stage 0 entry point — pMHC target PDB preparation.

Implements `specs/stage0_target_specification.md`. Cleans the primary
target (PDB 3HPJ, WT1 / HLA-A*02:01) and the positive-control target
(PDB 2BNR, Jenkins NY1-B04 / HLA-A*02:01) into canonical
chain-renamed, hetatm-stripped PDB files and emits a single
`target.yaml` manifest used by stage 1 (RFdiffusion).

Chain ordering convention (HC=A, β2m=B, peptide=C, binder=D appended in
stage 1) — HADRUP_JENKINS_2025. The binder chain is appended last with
a >12 Å residue-numbering gap to avoid AF2 chain-merge artefacts —
BAKER_LAB_2025 (stage 1's responsibility).
"""

from __future__ import annotations

import argparse
import shutil
import sys
from collections.abc import Iterable
from pathlib import Path
from typing import Literal

import yaml
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from pydantic import BaseModel, Field

THREE_TO_ONE: dict[str, str] = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

HC_LENGTH_RANGE = (260, 290)
B2M_LENGTH_RANGE = (94, 105)


class ChainMap(BaseModel):
    heavy_chain: str
    beta2m: str
    peptide: str
    binder: str


class BinderRange(BaseModel):
    length_min: int
    length_max: int


class TargetConfig(BaseModel):
    target_id: str
    pdb_id: str
    peptide_sequence: str
    mhc_allele: str
    chains: ChainMap
    hotspots_peptide: list[str] = Field(default_factory=list)
    hotspots_mhc: list[str] = Field(default_factory=list)
    binder: BinderRange
    cleaned_pdb_path: str


class ResolvedHotspot(BaseModel):
    token: str
    chain: str
    resnum_clean: int
    resnum_rcsb: int


class ResolvedChain(BaseModel):
    role: Literal["heavy_chain", "beta2m", "peptide", "binder"]
    length: int | None = None


class ResolvedTarget(BaseModel):
    target_id: str
    pdb_id: str
    cleaned_pdb: str
    peptide_sequence: str
    chains: dict[str, ResolvedChain]
    hotspots: dict[str, list[ResolvedHotspot]]
    binder: BinderRange


class TargetManifest(BaseModel):
    cycle: str
    mock: bool
    primary: ResolvedTarget
    positive_control: ResolvedTarget


def parse_hotspot(token: str) -> tuple[str, int]:
    if len(token) < 2 or not token[0].isalpha():
        raise ValueError(f"hotspot {token!r} must be <chain><resnum>, e.g. C1 or A155")
    chain_id = token[0]
    try:
        resnum = int(token[1:])
    except ValueError as exc:
        raise ValueError(f"hotspot {token!r} has non-integer resnum") from exc
    return chain_id, resnum


def load_target_config(path: Path) -> TargetConfig:
    with path.open() as fh:
        return TargetConfig.model_validate(yaml.safe_load(fh))


class _Keep(Select):
    """Reject HETATMs, waters, hydrogens, B-altlocs, and chains outside keep_set."""

    def __init__(self, keep_chains: Iterable[str]) -> None:
        self._keep = set(keep_chains)

    def accept_chain(self, chain: Chain) -> bool:
        return chain.id in self._keep

    def accept_residue(self, residue: Residue) -> bool:
        hetflag, _, _ = residue.id
        if hetflag.strip() != "":
            return False
        return residue.get_resname() not in {"HOH", "WAT"}

    def accept_atom(self, atom) -> bool:  # type: ignore[no-untyped-def]
        if atom.element == "H" or atom.get_name().startswith("H"):
            return False
        altloc = atom.get_altloc()
        return altloc in {" ", "", "A"}


def _residues(chain: Chain) -> list[Residue]:
    return [r for r in chain.get_residues() if r.id[0].strip() == ""]


def _rename_chains_safely(model: Model, mapping: dict[str, str]) -> None:
    """Rename chains via a two-pass temp-suffix to avoid Biopython's
    silent collision rename (spec pitfall).
    """
    for chain in list(model.get_chains()):
        if chain.id in mapping:
            chain.id = f"__tmp_{chain.id}"
    for chain in list(model.get_chains()):
        if chain.id.startswith("__tmp_"):
            original = chain.id.removeprefix("__tmp_")
            chain.id = mapping[original]


def _renumber_peptide(chain: Chain) -> dict[int, int]:
    """Renumber a peptide chain 1..N. Returns clean→rcsb resnum map."""
    residues = _residues(chain)
    clean_to_rcsb: dict[int, int] = {}
    for new_resnum, residue in enumerate(residues, start=1):
        clean_to_rcsb[new_resnum] = residue.id[1]
        residue.id = (" ", new_resnum, " ")
    return clean_to_rcsb


def _chain_one_letter(chain: Chain) -> str:
    return "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in _residues(chain))


def _prepend_remarks(pdb_path: Path, lines: list[str]) -> None:
    body = pdb_path.read_text()
    pdb_path.write_text("\n".join(lines) + "\n" + body)


def clean_target(
    raw_pdb: Path,
    cleaned_pdb_out: Path,
    cfg: TargetConfig,
) -> ResolvedTarget:
    """Clean one raw PDB into the canonical layout and return its manifest entry."""
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(cfg.pdb_id, str(raw_pdb))
    model: Model = structure[0]

    keep_originals = [cfg.chains.heavy_chain, cfg.chains.beta2m, cfg.chains.peptide]
    missing = [c for c in keep_originals if c not in {ch.id for ch in model.get_chains()}]
    if missing:
        raise ValueError(f"{cfg.pdb_id}: chains {missing} absent from raw PDB")

    rename_map = {
        cfg.chains.heavy_chain: "A",
        cfg.chains.beta2m: "B",
        cfg.chains.peptide: "C",
    }
    _rename_chains_safely(model, rename_map)

    hc_len = len(_residues(model["A"]))
    b2m_len = len(_residues(model["B"]))
    pep_len = len(_residues(model["C"]))
    if not HC_LENGTH_RANGE[0] <= hc_len <= HC_LENGTH_RANGE[1]:
        raise ValueError(f"{cfg.pdb_id}: HC chain A length {hc_len} outside {HC_LENGTH_RANGE}")
    if not B2M_LENGTH_RANGE[0] <= b2m_len <= B2M_LENGTH_RANGE[1]:
        raise ValueError(f"{cfg.pdb_id}: β2m chain B length {b2m_len} outside {B2M_LENGTH_RANGE}")
    if pep_len != len(cfg.peptide_sequence):
        raise ValueError(
            f"{cfg.pdb_id}: peptide chain C length {pep_len} != "
            f"config peptide_sequence length {len(cfg.peptide_sequence)}"
        )

    pep_clean_to_rcsb = _renumber_peptide(model["C"])

    pep_seq = _chain_one_letter(model["C"])
    if pep_seq != cfg.peptide_sequence:
        raise ValueError(
            f"{cfg.pdb_id}: cleaned peptide sequence {pep_seq!r} != "
            f"config {cfg.peptide_sequence!r}"
        )

    cleaned_pdb_out.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(cleaned_pdb_out), select=_Keep({"A", "B", "C"}))

    rcsb_chains = ", ".join(f"{src}->{dst}" for src, dst in rename_map.items())
    pep_orig = ", ".join(f"{c}->{r}" for c, r in pep_clean_to_rcsb.items())
    _prepend_remarks(
        cleaned_pdb_out,
        [
            f"REMARK 999 SOURCE_PDB {cfg.pdb_id}",
            f"REMARK 999 CHAIN_RENAME {rcsb_chains}",
            f"REMARK 999 PEPTIDE_ORIGINAL_NUMBERING {pep_orig}",
        ],
    )

    reparsed = PDBParser(QUIET=True).get_structure(cfg.pdb_id, str(cleaned_pdb_out))
    final_chains = {c.id for c in reparsed[0].get_chains()}
    if final_chains != {"A", "B", "C"}:
        raise ValueError(f"{cfg.pdb_id}: cleaned PDB has chains {final_chains}, expected {{A,B,C}}")

    resolved_hotspots = {
        "peptide": _resolve_hotspots(cfg.hotspots_peptide, reparsed[0], pep_clean_to_rcsb),
        "mhc": _resolve_hotspots(cfg.hotspots_mhc, reparsed[0], pep_clean_to_rcsb={}),
    }

    return ResolvedTarget(
        target_id=cfg.target_id,
        pdb_id=cfg.pdb_id,
        cleaned_pdb=str(cleaned_pdb_out),
        peptide_sequence=pep_seq,
        chains={
            "A": ResolvedChain(role="heavy_chain", length=hc_len),
            "B": ResolvedChain(role="beta2m", length=b2m_len),
            "C": ResolvedChain(role="peptide", length=pep_len),
            "D": ResolvedChain(role="binder"),
        },
        hotspots=resolved_hotspots,
        binder=cfg.binder,
    )


def _resolve_hotspots(
    tokens: list[str],
    model: Model,
    pep_clean_to_rcsb: dict[int, int],
) -> list[ResolvedHotspot]:
    resolved: list[ResolvedHotspot] = []
    for token in tokens:
        chain_id, resnum_clean = parse_hotspot(token)
        if chain_id not in model:
            raise ValueError(f"hotspot {token}: chain {chain_id} absent from cleaned PDB")
        chain = model[chain_id]
        if (" ", resnum_clean, " ") not in chain:
            raise ValueError(
                f"hotspot {token}: residue {resnum_clean} absent from chain {chain_id}"
            )
        rcsb = (
            pep_clean_to_rcsb.get(resnum_clean, resnum_clean) if chain_id == "C" else resnum_clean
        )
        resolved.append(
            ResolvedHotspot(
                token=token, chain=chain_id, resnum_clean=resnum_clean, resnum_rcsb=rcsb
            )
        )
    return resolved


def _mock_resolved(cfg: TargetConfig, cleaned_pdb_out: Path) -> ResolvedTarget:
    """Mock-mode ResolvedTarget — echoes config, no PDB inspection."""
    pep = [
        ResolvedHotspot(
            token=t,
            chain=parse_hotspot(t)[0],
            resnum_clean=parse_hotspot(t)[1],
            resnum_rcsb=parse_hotspot(t)[1],
        )
        for t in cfg.hotspots_peptide
    ]
    mhc = [
        ResolvedHotspot(
            token=t,
            chain=parse_hotspot(t)[0],
            resnum_clean=parse_hotspot(t)[1],
            resnum_rcsb=parse_hotspot(t)[1],
        )
        for t in cfg.hotspots_mhc
    ]
    return ResolvedTarget(
        target_id=cfg.target_id,
        pdb_id=cfg.pdb_id,
        cleaned_pdb=str(cleaned_pdb_out),
        peptide_sequence=cfg.peptide_sequence,
        chains={
            "A": ResolvedChain(role="heavy_chain"),
            "B": ResolvedChain(role="beta2m"),
            "C": ResolvedChain(role="peptide", length=len(cfg.peptide_sequence)),
            "D": ResolvedChain(role="binder"),
        },
        hotspots={"peptide": pep, "mhc": mhc},
        binder=cfg.binder,
    )


def _cycle_from_out(out: Path) -> str:
    for part in out.parts:
        if part.startswith("cycle_"):
            return part.removeprefix("cycle_")
    return "01"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary-config", required=True, type=Path)
    parser.add_argument("--control-config", required=True, type=Path)
    parser.add_argument("--primary-pdb-in", type=Path)
    parser.add_argument("--control-pdb-in", type=Path)
    parser.add_argument("--primary-pdb-out", required=True, type=Path)
    parser.add_argument("--control-pdb-out", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--mock", action="store_true")
    args = parser.parse_args(argv)

    primary_cfg = load_target_config(args.primary_config)
    control_cfg = load_target_config(args.control_config)
    cycle = _cycle_from_out(args.out)

    if args.mock:
        repo_root = Path(__file__).resolve().parents[2]
        fixtures = repo_root / "tests" / "fixtures"
        args.primary_pdb_out.parent.mkdir(parents=True, exist_ok=True)
        args.control_pdb_out.parent.mkdir(parents=True, exist_ok=True)
        if not args.primary_pdb_out.exists():
            shutil.copy(fixtures / "target_3hpj_clean.pdb", args.primary_pdb_out)
        if not args.control_pdb_out.exists():
            shutil.copy(fixtures / "target_2bnr_clean.pdb", args.control_pdb_out)
        manifest = TargetManifest(
            cycle=cycle,
            mock=True,
            primary=_mock_resolved(primary_cfg, args.primary_pdb_out),
            positive_control=_mock_resolved(control_cfg, args.control_pdb_out),
        )
    else:
        if args.primary_pdb_in is None or args.control_pdb_in is None:
            parser.error("--primary-pdb-in and --control-pdb-in required in non-mock mode")
        manifest = TargetManifest(
            cycle=cycle,
            mock=False,
            primary=clean_target(args.primary_pdb_in, args.primary_pdb_out, primary_cfg),
            positive_control=clean_target(args.control_pdb_in, args.control_pdb_out, control_cfg),
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        yaml.safe_dump(manifest.model_dump(), fh, sort_keys=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
