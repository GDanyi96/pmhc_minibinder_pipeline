"""Stage 6 entry point — render report.

See `specs/stage6_reporting.md` for the implementation contract.
This file is a placeholder; real logic lands in the stage-6 implementation session.
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
        print(f"[mock] render_report: config={args.config} out={args.out}")
        return 0
    raise NotImplementedError("stage 6 not implemented; see specs/stage6_reporting.md")


if __name__ == "__main__":
    sys.exit(main())
