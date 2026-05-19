# mypy: disable-error-code="no-untyped-call,attr-defined,no-any-return,misc"
"""DEPRECATED -- real-mode bypassed in cycle 02+.

Stage 2 designs (``scripts/run_stage2.py``) now invokes upstream
``dauparas/ProteinMPNN`` directly via the canonical multi-PDB batch
pattern (parse_multiple_chains -> assign_fixed_chains -> protein_mpnn_run),
mirroring ``scripts/run_stage1.py``'s direct-subprocess style. This
wrapper's ``run_real()`` carries three confirmed bugs:

* ``--fixed_chains`` / ``--designed_chain`` are not valid
  ``protein_mpnn_run.py`` flags (upstream uses ``--chain_id_jsonl``).
* ``assert_no_chain_d`` is misplaced relative to the splice step (binder
  IS chain D after splicing).
* Never validated against real data.

The mock paths remain functional and are still used by the cycle-1
controls regression. Do not invoke ``run_real`` from new code; do not
extend it. See ``specs/stage2_proteinmpnn_af2.md`` and the deprecation
notes in ``PIPELINE_STATUS.md`` for the rationale.

Mock mode: copies ``tests/fixtures/proteinmpnn/sample.fasta`` to
``<out_folder>/<pdb_stem>.fasta``. No GPU.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from Bio.PDB import PDBParser

PROTEINMPNN_DIR = Path(os.environ.get("PROTEINMPNN_DIR", "/workspace/ProteinMPNN"))
SAMPLE_FIXTURE = Path("tests/fixtures/proteinmpnn/sample.fasta")
DEFAULT_TEMP = 0.1
DEFAULT_NUM_SEQ = 10
DEFAULT_SEED = 42


def assert_no_chain_d(pdb_path: Path) -> None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    chains = [c.id for c in structure[0].get_chains()]
    if "D" in chains:
        raise AssertionError(
            f"{pdb_path}: chain D already present " "(expected fixed A+B+C scaffold from Stage 0)"
        )


def run_mock(pdb_path: Path, out_folder: Path) -> None:
    out_folder.mkdir(parents=True, exist_ok=True)
    dest = out_folder / f"{pdb_path.stem}.fasta"
    shutil.copy(SAMPLE_FIXTURE, dest)


def run_real(
    pdb_path: Path,
    chain_id_jsonl: Path,
    fixed_chains: list[str],
    designed_chain: str,
    sampling_temp: float,
    num_seq_per_target: int,
    out_folder: Path,
    seed: int,
) -> None:
    out_folder.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(PROTEINMPNN_DIR / "protein_mpnn_run.py"),
        "--pdb_path",
        str(pdb_path),
        "--chain_id_jsonl",
        str(chain_id_jsonl),
        "--fixed_chains",
        *fixed_chains,
        "--designed_chain",
        designed_chain,
        "--sampling_temp",
        str(sampling_temp),
        "--num_seq_per_target",
        str(num_seq_per_target),
        "--out_folder",
        str(out_folder),
        "--seed",
        str(seed),
    ]
    subprocess.run(cmd, check=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb_path", type=Path)
    parser.add_argument("--chain_id_jsonl", type=Path)
    parser.add_argument("--fixed_chains", nargs="+", default=["A", "B", "C"])
    parser.add_argument("--designed_chain", default="D")
    parser.add_argument("--sampling_temp", type=float, default=DEFAULT_TEMP)
    parser.add_argument("--num_seq_per_target", type=int, default=DEFAULT_NUM_SEQ)
    parser.add_argument("--out_folder", type=Path)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--mock", action="store_true")
    # Back-compat with the original stub signature used by the existing rule.
    parser.add_argument("--config")
    parser.add_argument("--out")
    args = parser.parse_args(argv)

    pdb_path: Path | None = args.pdb_path or (Path(args.config) if args.config else None)
    out_folder: Path | None = args.out_folder or (Path(args.out).parent if args.out else None)
    if pdb_path is None or out_folder is None:
        parser.error("must provide --pdb_path and --out_folder (or legacy --config/--out)")

    if args.mock:
        run_mock(pdb_path, out_folder)
        print(f"[mock] run_proteinmpnn: copied sample fasta to {out_folder}")
        return 0

    assert_no_chain_d(pdb_path)
    if args.chain_id_jsonl is None:
        parser.error("--chain_id_jsonl required in real mode")
    run_real(
        pdb_path,
        args.chain_id_jsonl,
        args.fixed_chains,
        args.designed_chain,
        args.sampling_temp,
        args.num_seq_per_target,
        out_folder,
        args.seed,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
