"""Populate N1 (scrambled P1) and N2 (random natural-frequency) binder sequences.

Seeded with `random.Random(42)`. Idempotent: only writes entries whose
`binder_sequence` is currently null. The manifest is round-tripped through
PyYAML with `sort_keys=False` to preserve field order.

This script is deterministic and has no GPU/Docker dependency; `--mock`
behaves identically to real mode and is provided only for CLAUDE.md uniformity.
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path
from typing import Any

import yaml

DEFAULT_MANIFEST = Path("data/controls/controls_manifest.yaml")
SEED = 42


def _load_manifest(path: Path) -> dict[str, Any]:
    with path.open() as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict) or "controls" not in data:
        raise ValueError(f"{path}: missing top-level 'controls' list")
    return data


def _dump_manifest(path: Path, data: dict[str, Any]) -> None:
    with path.open("w") as fh:
        yaml.safe_dump(data, fh, sort_keys=False, default_flow_style=False, width=1000)


def _control(data: dict[str, Any], control_id: str) -> dict[str, Any]:
    for entry in data["controls"]:
        if entry.get("id") == control_id:
            if not isinstance(entry, dict):
                raise TypeError(f"control {control_id}: entry is not a mapping")
            return entry
    raise KeyError(f"control {control_id!r} not found in manifest")


def scramble(sequence: str, rng: random.Random) -> str:
    chars = list(sequence)
    rng.shuffle(chars)
    return "".join(chars)


def random_sequence(length: int, frequencies: dict[str, float], rng: random.Random) -> str:
    alphabet = sorted(frequencies)
    weights = [frequencies[a] for a in alphabet]
    return "".join(rng.choices(alphabet, weights=weights, k=length))


def populate(data: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    """Fill any null N1/N2 binder_sequence in-place. Returns (data, written_ids)."""
    p1 = _control(data, "P1")
    p1_seq = p1.get("binder_sequence")
    if not isinstance(p1_seq, str):
        raise ValueError("P1 binder_sequence must be a string before negatives can be generated")

    rng = random.Random(SEED)
    written: list[str] = []

    n1 = _control(data, "N1")
    if n1.get("binder_sequence") is None:
        n1_length = int(n1["binder_length"])
        if n1_length != len(p1_seq):
            raise ValueError(f"N1 binder_length {n1_length} != P1 sequence length {len(p1_seq)}")
        n1["binder_sequence"] = scramble(p1_seq, rng)
        written.append("N1")

    n2 = _control(data, "N2")
    if n2.get("binder_sequence") is None:
        n2_length = int(n2["binder_length"])
        frequencies = n2.get("aa_frequencies")
        if not isinstance(frequencies, dict):
            raise ValueError("N2 aa_frequencies block missing or malformed")
        freqs_typed: dict[str, float] = {str(k): float(v) for k, v in frequencies.items()}
        n2["binder_sequence"] = random_sequence(n2_length, freqs_typed, rng)
        written.append("N2")

    return data, written


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--mock", action="store_true")
    args = parser.parse_args(argv)

    data = _load_manifest(args.config)
    data, written = populate(data)
    if written:
        _dump_manifest(args.config, data)
        print(f"generate_negatives: wrote {written} into {args.config}")
    else:
        print(f"generate_negatives: N1 and N2 already populated in {args.config}; no-op")
    if args.mock:
        print("[mock] generate_negatives runs identically in mock and real modes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
