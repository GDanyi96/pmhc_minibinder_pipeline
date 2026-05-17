"""Stage 4 entry point — embed designs.

See `specs/stage4_diversity.md` for the implementation contract.
This file is a placeholder; real logic lands in the stage-4 implementation session.
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
        print(f"[mock] embed_designs: config={args.config} out={args.out}")
        return 0
    raise NotImplementedError("stage 4 not implemented; see specs/stage4_diversity.md")


if __name__ == "__main__":
    sys.exit(main())
