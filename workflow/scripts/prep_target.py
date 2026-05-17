"""Stage 0 entry point — target PDB preparation.

See `specs/stage0_target_specification.md` for the implementation contract.
This file is a placeholder; real logic lands in the stage-0 implementation
session.
"""

from __future__ import annotations

import argparse
import sys


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--mock", action="store_true")
    args = parser.parse_args(argv)
    if args.mock:
        print(f"[mock] prep_target: config={args.config} out={args.out}")
        return 0
    raise NotImplementedError("stage 0 not implemented; see specs/stage0_target_specification.md")


if __name__ == "__main__":
    sys.exit(main())
