"""Stage 3 entry point — crosspan.

See `specs/stage3_crosspan.md` for the implementation contract.
This file is a placeholder; real logic lands in the stage-3 implementation session.
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
        print(f"[mock] crosspan: config={args.config} out={args.out}")
        return 0
    raise NotImplementedError("stage 3 not implemented; see specs/stage3_crosspan.md")


if __name__ == "__main__":
    sys.exit(main())
