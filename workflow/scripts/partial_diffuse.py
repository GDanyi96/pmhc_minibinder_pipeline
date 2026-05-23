# mypy: disable-error-code="no-untyped-call,attr-defined,no-any-return,misc"
"""Cycle 03 -- RFdiffusion partial-diffusion CLI wrapper.

Partial diffusion noises an existing backbone (a BAKER scaffold for sub-run
A, our cycle-02 hero for sub-run B) by ``partial_T`` reverse steps and
re-denoises it, remodeling the binder while the target chains stay fixed as
the motif. This is BAKER_LAB_2025's privileged-scaffold-reuse pathway, and
the cycle-03 answer to cycle-02's placement deficit (Trap #30 context).

Canonical args mirror RFdiffusion's partial-diffusion examples
(``/workspace/RFdiffusion/examples/partial_*``):
``diffuser.partial_T=<N>``, ``inference.input_pdb=<seed>``,
``denoiser.noise_scale_ca=0`` for binder mode. ``partial_T`` and
``noise_scale_ca`` are read from configs/rfdiffusion_subrun_{a,b}.yaml;
defaults ``partial_T=10, noise_scale_ca=0`` are a conservative blind default
(no published optimum for pMHC scaffold reuse -- recalibrate in cycle 04).

Real mode shells out to ``run_inference.py`` via the SE3nv interpreter
(see run_rfdiffusion._resolve_rfdiffusion_python) and redirects Hydra's run
dir to ``<out_dir>/hydra_logs``. Mock mode synthesizes a geometry-passing
4-chain design (A/B/C motif copied from the reference + a chain-D binder
backbone through the hotspot-cluster centroid) so the merge halt gate is
exercised end-to-end with no GPU.
"""

from __future__ import annotations

import os

os.environ.pop("LD_LIBRARY_PATH", None)

import argparse
import logging
import subprocess
import sys
from pathlib import Path

from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from workflow.scripts.run_rfdiffusion import (
    RFDIFFUSION_DIR,
    _resolve_rfdiffusion_python,
)

logger = logging.getLogger("partial_diffuse")
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

DEFAULT_PARTIAL_T = 10
DEFAULT_NOISE_SCALE_CA = 0
_MOCK_BINDER_CHAIN = "D"
_MOTIF_CHAINS = ("A", "B", "C")

# Sub-run A aligned scaffolds carry the target as the BAKER-format truncated
# motif: chain B = HLA[1:180], chain C = peptide[1:9] (align_baker_scaffolds).
# These two lengths are static (fixed by the truncated target geometry); only
# the chain-A binder length varies per scaffold and is derived per-design.
_TRUNCATED_HLA_LEN = 180
_TRUNCATED_PEPTIDE_LEN = 9


def _motif_residues(reference_pdb: Path) -> dict[str, list]:  # type: ignore[type-arg]
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(reference_pdb.stem, str(reference_pdb))
    model = structure[0]
    out: dict[str, list] = {}  # type: ignore[type-arg]
    for cid in _MOTIF_CHAINS:
        if cid not in model:
            raise ValueError(f"{reference_pdb}: missing motif chain {cid!r}")
        out[cid] = [r for r in model[cid].get_residues() if r.id[0].strip() == ""]
    return out


def _binder_backbone_chain(n: int) -> Chain:
    """Chain-D backbone through the hotspot-cluster centroid (mock).

    Same straight-backbone formula as the Stage 1 / cycle-03 fixtures so the
    synthesized design passes the geometry gate (>=3 CA within 10 A of a
    hotspot, length in range, no internal clash).
    """
    import numpy as np
    from Bio.PDB.Atom import Atom
    from Bio.PDB.Residue import Residue

    chain = Chain(_MOCK_BINDER_CHAIN)
    for i in range(n):
        x = 1.5
        y = i * 3.9 - (n - 1) * 3.9 / 2 + 3.5
        z = 2.5
        residue = Residue((" ", i + 1, " "), "ALA", "")
        coord = np.array([x, y, z], dtype=float)
        atom = Atom("CA", coord, 0.0, 1.0, " ", "CA", i + 1, "C")
        residue.add(atom)
        chain.add(residue)
    return chain


def _synthesize_mock_design(out_pdb: Path, reference_pdb: Path, binder_length: int) -> None:
    """Write a geometry-passing 4-chain design (A/B/C motif + D binder)."""
    motif = _motif_residues(reference_pdb)
    out_structure = Structure(out_pdb.stem)
    out_model = Model(0)
    out_structure.add(out_model)
    for cid in _MOTIF_CHAINS:
        chain = Chain(cid)
        for residue in motif[cid]:
            chain.add(residue.copy())
        out_model.add(chain)
    out_model.add(_binder_backbone_chain(binder_length))

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(out_structure)
    io.save(str(out_pdb))


def _derive_contigs_subrun_a(aligned_scaffold_pdb: Path) -> str:
    """Build the RFdiffusion contig string for one aligned sub-run A scaffold.

    The aligned scaffold (align_baker_scaffolds output) is laid out chain
    A = binder, B = HLA[1:180], C = peptide[1:9]. The binder length N is the
    chain-A standard-residue CA count.

    RFdiffusion requires ``contigmap.contigs`` in partial-diffusion mode too
    (it is the residue mask, not an optional bias -- Trap #32). The binder
    slot is a *bare* ``N-N`` range (unprefixed = redesigned), while ``B``/``C``
    are letter-prefixed motifs (preserved). A letter on the binder slot would
    tell RFdiffusion to preserve chain A, defeating partial diffusion. See
    RFdiffusion README "Partial diffusion": e.g. ``[100-100/0 B1-150]``.
    """
    parser = PDBParser(QUIET=True)
    structure: Structure = parser.get_structure(
        aligned_scaffold_pdb.stem, str(aligned_scaffold_pdb)
    )
    model = structure[0]
    if "A" in model:
        n = sum(1 for r in model["A"].get_residues() if r.id[0].strip() == "" and "CA" in r)
    else:
        n = 0
    if n == 0:
        raise ValueError(
            f"partial_diffuse: contigs not derivable from {aligned_scaffold_pdb} "
            "— chain A empty or malformed"
        )
    return f"[{n}-{n}/0 B1-{_TRUNCATED_HLA_LEN}/0 C1-{_TRUNCATED_PEPTIDE_LEN}]"


def partial_diffuse_one(
    input_pdb: Path,
    out_prefix: Path,
    seed: int,
    ckpt: Path,
    partial_T: int,
    noise_scale_ca: int,
    hydra_dir: Path,
) -> None:
    """Run one RFdiffusion partial-diffusion design (real mode)."""
    rfdiffusion_python = _resolve_rfdiffusion_python()
    hydra_dir.mkdir(parents=True, exist_ok=True)
    contig = _derive_contigs_subrun_a(input_pdb)
    cmd = [
        rfdiffusion_python,
        str(RFDIFFUSION_DIR / "scripts" / "run_inference.py"),
        f"inference.input_pdb={input_pdb}",
        f"inference.output_prefix={out_prefix}",
        "inference.num_designs=1",
        f"inference.ckpt_override_path={ckpt}",
        f"diffuser.partial_T={partial_T}",
        f"denoiser.noise_scale_ca={noise_scale_ca}",
        "denoiser.noise_scale_frame=0",
        "inference.deterministic=True",
        f"inference.design_startnum={seed}",
        f"hydra.run.dir={hydra_dir}",
        f"contigmap.contigs={contig}",
    ]
    logger.info("RFdiffusion partial diffusion: %s", cmd)
    subprocess.run(cmd, check=True, env={**os.environ})


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-pdb", type=Path, required=True)
    parser.add_argument("--out-pdb", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--partial-t", type=int, default=DEFAULT_PARTIAL_T)
    parser.add_argument("--noise-scale-ca", type=int, default=DEFAULT_NOISE_SCALE_CA)
    parser.add_argument("--ckpt", type=Path, default=None)
    parser.add_argument("--reference-pdb", type=Path, default=None)
    parser.add_argument("--binder-length", type=int, default=80)
    parser.add_argument("--hydra-dir", type=Path, default=None)
    parser.add_argument("--mock", action="store_true")
    args = parser.parse_args(argv)

    if args.mock:
        reference = args.reference_pdb or Path("tests/fixtures/stage1/mock_clean.pdb")
        _synthesize_mock_design(args.out_pdb, reference, args.binder_length)
        logger.info("mock: synthesized %s", args.out_pdb)
        return 0

    if args.ckpt is None:
        parser.error("--ckpt is required in real mode")
    partial_diffuse_one(
        input_pdb=args.input_pdb,
        out_prefix=args.out_pdb.with_suffix(""),
        seed=args.seed,
        ckpt=args.ckpt,
        partial_T=args.partial_t,
        noise_scale_ca=args.noise_scale_ca,
        hydra_dir=args.hydra_dir or args.out_pdb.parent / "hydra_logs",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
