"""Stage 2b — ColabFold (AF2-multimer) wrapper.

Real mode: shells out to the Docker image `colabfold/colabfold:1.5.5-cuda12.2.2`
with `--num-recycle 6 --model-type alphafold2_multimer_v3 --use-gpu-relax`
per MARES_IOANNIDIS_2025.

Mock mode: copies pre-baked fixtures from `tests/fixtures/stage2/controls_colabfold/`
into per-id subdirectories of `--out-dir`. No Docker, no GPU.

FASTA contract: one entry per id, header `>{id}`, body
`<HC>:<beta2m>:<peptide>:<binder>` (colon-separated multimer). Chain order
A=HC, B=beta2m, C=peptide, D=binder — pitfall #1 in specs/stage2_proteinmpnn_af2.md.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

DOCKER_IMAGE = "colabfold/colabfold:1.5.5-cuda12.2.2"
DEFAULT_NUM_RECYCLES = 6
FIXTURES_DIR = Path("tests/fixtures/stage2/controls_colabfold")


def _ids_from_fasta_dir(fasta_dir: Path) -> list[str]:
    ids: list[str] = []
    for fasta in sorted(fasta_dir.glob("*.fasta")):
        for line in fasta.read_text().splitlines():
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    if not ids:
        raise ValueError(f"no FASTA entries found under {fasta_dir}")
    return ids


def run_mock(fasta_dir: Path, out_dir: Path, fixtures: Path = FIXTURES_DIR) -> None:
    ids = _ids_from_fasta_dir(fasta_dir)
    for control_id in ids:
        dest = out_dir / control_id
        dest.mkdir(parents=True, exist_ok=True)
        copied = 0
        for src in fixtures.glob(f"{control_id}_*_rank_001*"):
            shutil.copy(src, dest / src.name)
            copied += 1
        if copied < 2:
            raise FileNotFoundError(
                f"mock fixtures missing for {control_id}: "
                f"expected scores.json + pdb under {fixtures}"
            )


def run_real(fasta_dir: Path, out_dir: Path, num_recycles: int) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "docker",
        "run",
        "--gpus",
        "all",
        "-v",
        "/workspace:/workspace",
        DOCKER_IMAGE,
        "colabfold_batch",
        "--num-recycle",
        str(num_recycles),
        "--model-type",
        "alphafold2_multimer_v3",
        "--use-gpu-relax",
        str(fasta_dir),
        str(out_dir),
    ]
    subprocess.run(cmd, check=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta-dir", type=Path)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--num-recycles", type=int, default=DEFAULT_NUM_RECYCLES)
    parser.add_argument("--fixtures-dir", type=Path, default=FIXTURES_DIR)
    parser.add_argument("--mock", action="store_true")
    # Back-compat with the original stub signature used by the existing
    # Snakemake rule. Removed once Part B rewires rule 03.
    parser.add_argument("--config")
    parser.add_argument("--out")
    args = parser.parse_args(argv)

    fasta_dir: Path | None = args.fasta_dir
    if fasta_dir is None and args.config:
        fasta_dir = Path(args.config).parent
    out_dir: Path | None = args.out_dir
    if out_dir is None and args.out:
        out_dir = Path(args.out).parent
    if fasta_dir is None or out_dir is None:
        parser.error("must provide --fasta-dir and --out-dir (or legacy --config/--out)")

    if args.mock:
        run_mock(fasta_dir, out_dir, fixtures=args.fixtures_dir)
        print(f"[mock] run_colabfold: copied fixtures for {len(list(out_dir.iterdir()))} ids")
        return 0

    run_real(fasta_dir, out_dir, args.num_recycles)
    return 0


if __name__ == "__main__":
    sys.exit(main())
