# mypy: disable-error-code="no-untyped-call,attr-defined,no-any-return,misc"
"""Pod-side one-shot setup for Cycle 03 real inputs. NOT run by CI / mock.

Cycle 03 depends on two artefacts that live only on the RunPod pod and are
gitignored, so they are materialized here after the PR is merged and pulled:

    python workflow/scripts/setup_cycle03_inputs.py

1. Sub-run A inputs -- symlink BAKER's published scaffold library into the
   repo so configs/rfdiffusion_subrun_a.yaml's glob resolves:
       /workspace/pMHCI_binder_design/scaffolds/scaf*.pdb
         -> data/scaffolds/baker_library/scaf*.pdb

2. Sub-run B seed -- the cycle-02 hero binder (design_2079_seq00, chain D of
   the AF2 prediction) stitched onto the canonical 3hpj_clean.pdb as a full
   4-chain complex (A=HC, B=beta2m, C=peptide, D=binder). RFdiffusion partial
   diffusion needs the target context as the fixed motif, so we keep the
   whole complex, not the binder alone. Reuses
   workflow.scripts.splice_binder (which takes chain D from the prediction
   and A/B/C from the canonical reference), guaranteeing the same chain frame
   as sub-run A's aligned scaffolds:
       results/cycle_02/.../design_2079_seq00_unrelaxed_rank_001_*.pdb (chain D)
       + data/targets/3hpj_clean.pdb (chains A/B/C)
         -> data/seeds/design_2079_binder.pdb

Idempotent: existing outputs are left untouched.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from workflow.scripts import splice_binder

logger = logging.getLogger("setup_cycle03_inputs")
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

DEFAULT_BAKER_SCAFFOLDS = Path("/workspace/pMHCI_binder_design/scaffolds")
DEFAULT_HERO_PRED_GLOB = (
    "results/cycle_02/stage2/af2_designs/predictions/design_2079_seq00/"
    "design_2079_seq00_unrelaxed_rank_001_*.pdb"
)
DEFAULT_REFERENCE_PDB = Path("data/targets/3hpj_clean.pdb")
DEFAULT_SCAFFOLD_DEST = Path("data/scaffolds/baker_library")
DEFAULT_SEED_OUT = Path("data/seeds/design_2079_binder.pdb")


def symlink_scaffolds(src_dir: Path, dest_dir: Path) -> int:
    if not src_dir.exists():
        raise FileNotFoundError(f"BAKER scaffold dir not found: {src_dir}")
    dest_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for scaffold in sorted(src_dir.glob("scaf*.pdb")):
        link = dest_dir / scaffold.name
        if link.exists() or link.is_symlink():
            continue
        link.symlink_to(scaffold.resolve())
        n += 1
    logger.info("symlinked %d scaffolds into %s", n, dest_dir)
    return n


def build_seed_complex(hero_glob: str, reference_pdb: Path, seed_out: Path) -> None:
    if seed_out.exists():
        logger.info("seed already present, skipping: %s", seed_out)
        return
    base = Path(hero_glob)
    matches = sorted(base.parent.glob(base.name))
    if not matches:
        raise FileNotFoundError(f"hero prediction not found: {hero_glob}")
    if not reference_pdb.exists():
        raise FileNotFoundError(f"reference {reference_pdb} missing; run prep_target.py first")
    seed_out.parent.mkdir(parents=True, exist_ok=True)
    splice_binder.splice_binder_onto_pmhc(matches[0], reference_pdb, seed_out)
    logger.info("wrote stitched seed complex: %s", seed_out)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baker-scaffolds", type=Path, default=DEFAULT_BAKER_SCAFFOLDS)
    parser.add_argument("--scaffold-dest", type=Path, default=DEFAULT_SCAFFOLD_DEST)
    parser.add_argument("--hero-glob", default=DEFAULT_HERO_PRED_GLOB)
    parser.add_argument("--reference-pdb", type=Path, default=DEFAULT_REFERENCE_PDB)
    parser.add_argument("--seed-out", type=Path, default=DEFAULT_SEED_OUT)
    args = parser.parse_args(argv)

    symlink_scaffolds(args.baker_scaffolds, args.scaffold_dest)
    build_seed_complex(args.hero_glob, args.reference_pdb, args.seed_out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
