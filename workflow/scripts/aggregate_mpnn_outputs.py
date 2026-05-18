"""Stage 2 designs — aggregate per-design ProteinMPNN FASTA outputs.

The ProteinMPNN driver (`workflow.scripts.run_proteinmpnn`) writes one
FASTA per backbone PDB (the upstream convention). For the Stage 2
designs fan-in we need a single ranked-by-NLL records file across all
200 backbones x 4 sequences = 800 records. This module walks the
per-design output directories, parses each FASTA header for the score
field, and emits a unified `sequences.jsonl`.

Real ProteinMPNN headers look like
`>T=0.1, sample=1, score=0.8543, global_score=0.852, ...`. The mock
fixture uses the simpler `>design_NNNNN_seqK mock-mpnn nll=1.234`. Both
are supported.

Seed helpers `_compute_mpnn_seed` and `_compute_af2_seed` re-derive
seeds from the formulas in `configs/seeds.yaml` and assert membership in
the corresponding per-stage reserved range — same pattern as
`workflow.scripts.run_rfdiffusion._compute_seed`.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import yaml
from pydantic import BaseModel

_MPNN_FORMULA = "cycle * 1000 + 200 + design_index * 4 + seq_index"
_AF2_FORMULA = "cycle * 1000 + 1000 + fan_in_rank"

_FLOAT = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
_SCORE_RE = re.compile(rf"score=({_FLOAT})")
_NLL_RE = re.compile(rf"nll=({_FLOAT})")


class _Seeds(BaseModel):
    formulas: dict[str, str]
    reserved: dict[str, list[int]]


def _load_seeds(seeds_yaml: Path) -> _Seeds:
    with seeds_yaml.open() as fh:
        return _Seeds.model_validate(yaml.safe_load(fh))


def _compute_mpnn_seed(cycle: int, design_index: int, seq_index: int, seeds: _Seeds) -> int:
    actual = seeds.formulas.get("proteinmpnn")
    if actual != _MPNN_FORMULA:
        raise ValueError(
            f"unsupported proteinmpnn seed formula {actual!r}; " f"expected {_MPNN_FORMULA!r}"
        )
    seed = cycle * 1000 + 200 + design_index * 4 + seq_index
    key = f"cycle_{cycle:02d}_proteinmpnn"
    if key not in seeds.reserved:
        raise ValueError(f"seeds.yaml: no reserved range for {key}")
    lo, hi = seeds.reserved[key]
    if not lo <= seed <= hi:
        raise ValueError(
            f"mpnn seed {seed} (cycle={cycle}, design_index={design_index}, "
            f"seq_index={seq_index}) outside reserved range [{lo}, {hi}]"
        )
    return seed


def _compute_af2_seed(cycle: int, fan_in_rank: int, seeds: _Seeds) -> int:
    actual = seeds.formulas.get("af2_fanin")
    if actual != _AF2_FORMULA:
        raise ValueError(
            f"unsupported af2_fanin seed formula {actual!r}; expected {_AF2_FORMULA!r}"
        )
    seed = cycle * 1000 + 1000 + fan_in_rank
    key = f"cycle_{cycle:02d}_af2_fanin"
    if key not in seeds.reserved:
        raise ValueError(f"seeds.yaml: no reserved range for {key}")
    lo, hi = seeds.reserved[key]
    if not lo <= seed <= hi:
        raise ValueError(
            f"af2 fan-in seed {seed} (cycle={cycle}, fan_in_rank={fan_in_rank}) "
            f"outside reserved range [{lo}, {hi}]"
        )
    return seed


def _parse_score(header: str) -> float:
    """Extract the per-residue NLL score from a ProteinMPNN FASTA header.

    Prefers the real-format `score=N` field; falls back to the mock-format
    `nll=N`. Raises if neither is present.
    """
    m = _SCORE_RE.search(header)
    if m is None:
        m = _NLL_RE.search(header)
    if m is None:
        raise ValueError(f"FASTA header missing score=/nll= field: {header!r}")
    return float(m.group(1))


def _parse_fasta(fasta_path: Path) -> list[tuple[str, str]]:
    """Return list of (header_line_without_leading_gt, sequence) pairs.

    ProteinMPNN's first record is the input sequence (the fixed scaffold);
    subsequent records are the sampled designs. We return all records and
    let the caller decide.
    """
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []
    for raw in fasta_path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts)))
            header = line[1:]
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def aggregate(
    per_design_fastas: list[Path],
    design_records: list[dict[str, str | int]],
    cycle: int,
    seeds: _Seeds,
    out_jsonl: Path,
) -> Path:
    """Walk per-design FASTAs, emit unified sequences.jsonl.

    `design_records[i]` corresponds to the i-th element of
    `per_design_fastas`; must carry at least `design_id`, `backbone_pdb`,
    `spliced_pdb`. ProteinMPNN's first FASTA record (the fixed scaffold)
    is skipped; only sampled designs go into the output.
    """
    if len(per_design_fastas) != len(design_records):
        raise ValueError(
            f"length mismatch: {len(per_design_fastas)} fastas vs "
            f"{len(design_records)} design records"
        )
    out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with out_jsonl.open("w") as fh:
        for design_index, (fasta_path, design_record) in enumerate(
            zip(per_design_fastas, design_records, strict=True)
        ):
            fasta_records = _parse_fasta(fasta_path)
            samples = fasta_records[1:]
            for seq_index, (header, seq) in enumerate(samples):
                nll = _parse_score(header)
                seed = _compute_mpnn_seed(cycle, design_index, seq_index, seeds)
                out_record = {
                    "design_id": str(design_record["design_id"]),
                    "seq_id": seq_index,
                    "seq": seq,
                    "mpnn_nll": nll,
                    "mpnn_seed": seed,
                    "backbone_pdb": str(design_record["backbone_pdb"]),
                    "spliced_pdb": str(design_record["spliced_pdb"]),
                }
                fh.write(json.dumps(out_record, sort_keys=True) + "\n")
    return out_jsonl
