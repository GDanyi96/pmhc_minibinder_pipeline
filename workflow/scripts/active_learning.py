"""Stage 5 entry point — active learning.

See `specs/stage5_active_learning.md` for the implementation contract.
This file is a placeholder; real logic lands in the stage-5 implementation session.
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
        print(f"[mock] active_learning: config={args.config} out={args.out}")
        return 0
    raise NotImplementedError("stage 5 not implemented; see specs/stage5_active_learning.md")


if __name__ == "__main__":
    sys.exit(main())
