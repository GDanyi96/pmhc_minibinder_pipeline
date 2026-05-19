"""Unit tests for workflow.scripts.aggregate_mpnn_outputs."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from workflow.scripts.aggregate_mpnn_outputs import (
    _AF2_FORMULA,
    _MPNN_FORMULA,
    _compute_af2_seed,
    _compute_mpnn_seed,
    _parse_score,
    _Seeds,
    aggregate,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
SEEDS_YAML = REPO_ROOT / "configs" / "seeds.yaml"


def _load_real_seeds() -> _Seeds:
    return _Seeds.model_validate(yaml.safe_load(SEEDS_YAML.read_text()))


def test_parse_score_real_format() -> None:
    header = "T=0.1, sample=1, score=0.8543, global_score=0.852, model_name=v_48_020"
    assert _parse_score(header) == pytest.approx(0.8543)


def test_parse_score_mock_format() -> None:
    header = "design_00000_seq0 mock-mpnn nll=1.234"
    assert _parse_score(header) == pytest.approx(1.234)


def test_parse_score_raises_when_no_field() -> None:
    with pytest.raises(ValueError, match="missing score=/nll="):
        _parse_score("design_X no_score_here")


def test_mpnn_seed_formula_matches_yaml() -> None:
    seeds = _load_real_seeds()
    assert seeds.formulas["proteinmpnn"] == _MPNN_FORMULA
    assert seeds.formulas["af2_fanin"] == _AF2_FORMULA


def test_mpnn_seed_cycle2_start() -> None:
    seeds = _load_real_seeds()
    assert _compute_mpnn_seed(2, 0, 0, seeds) == 2200


def test_mpnn_seed_cycle2_end() -> None:
    seeds = _load_real_seeds()
    assert _compute_mpnn_seed(2, 199, 3, seeds) == 2999


def test_af2_seed_cycle2_range() -> None:
    seeds = _load_real_seeds()
    assert _compute_af2_seed(2, 0, seeds) == 3000
    assert _compute_af2_seed(2, 99, seeds) == 3099


def test_seeds_unique_across_all_records() -> None:
    """All 800 MPNN records (200 designs * 4 seqs) must have unique seeds in
    the cycle 2 reserved range; all 100 AF2 fan-in records likewise."""
    seeds = _load_real_seeds()
    mpnn_seeds = {
        _compute_mpnn_seed(2, design_index, seq_index, seeds)
        for design_index in range(200)
        for seq_index in range(4)
    }
    assert len(mpnn_seeds) == 800
    mpnn_lo, mpnn_hi = seeds.reserved["cycle_02_proteinmpnn"]
    assert all(mpnn_lo <= s <= mpnn_hi for s in mpnn_seeds)

    af2_seeds = {_compute_af2_seed(2, rank, seeds) for rank in range(100)}
    assert len(af2_seeds) == 100
    af2_lo, af2_hi = seeds.reserved["cycle_02_af2_fanin"]
    assert all(af2_lo <= s <= af2_hi for s in af2_seeds)


def test_mpnn_seed_rejects_out_of_range() -> None:
    seeds = _load_real_seeds()
    with pytest.raises(ValueError, match="outside reserved range"):
        _compute_mpnn_seed(2, 200, 0, seeds)


def test_mpnn_seed_rejects_unknown_cycle() -> None:
    seeds = _load_real_seeds()
    with pytest.raises(ValueError, match="no reserved range"):
        _compute_mpnn_seed(7, 0, 0, seeds)


def test_aggregate_emits_records_sorted_by_caller(tmp_path: Path) -> None:
    """aggregate() emits records in input order; ranking is the caller's job."""
    fasta_a = tmp_path / "design_00000.fa"
    fasta_a.write_text(
        ">scaffold_fixed\nGGGG\n"
        ">T=0.1, sample=1, score=0.95\nMAEELA\n"
        ">T=0.1, sample=2, score=0.85\nMAEELB\n"
    )
    fasta_b = tmp_path / "design_00001.fa"
    fasta_b.write_text(">scaffold_fixed\nGGGG\n" ">T=0.1, sample=1, score=0.75\nMAEELC\n")
    out = tmp_path / "sequences.jsonl"
    seeds = _load_real_seeds()
    aggregate(
        per_design_fastas=[fasta_a, fasta_b],
        design_records=[
            {
                "design_id": "design_00000",
                "backbone_pdb": "stage1/design_00000.pdb",
                "spliced_pdb": "spliced/design_00000.pdb",
            },
            {
                "design_id": "design_00001",
                "backbone_pdb": "stage1/design_00001.pdb",
                "spliced_pdb": "spliced/design_00001.pdb",
            },
        ],
        cycle=2,
        seeds=seeds,
        out_jsonl=out,
    )
    records = [json.loads(line) for line in out.read_text().splitlines()]
    assert len(records) == 3
    assert records[0]["design_id"] == "design_00000" and records[0]["seq_id"] == 0
    assert records[1]["design_id"] == "design_00000" and records[1]["seq_id"] == 1
    assert records[2]["design_id"] == "design_00001" and records[2]["seq_id"] == 0
    assert records[0]["mpnn_nll"] == pytest.approx(0.95)
    assert records[0]["mpnn_seed"] == 2200
    assert records[1]["mpnn_seed"] == 2201
    assert records[2]["mpnn_seed"] == 2204


def test_aggregate_slices_binder_when_binder_length_provided(tmp_path: Path) -> None:
    """Real ProteinMPNN emits the full A+B+C+D complex sequence per
    record; the aggregator must keep only the last N residues when
    ``binder_length`` is supplied. Trap #6.
    """
    fasta = tmp_path / "design_00000.fa"
    fasta.write_text(
        ">scaffold_full\nAAABBBCCCDDD\n"
        ">T=0.1, sample=1, score=0.50, global_score=0.60\nAAABBBCCCEEE\n"
    )
    out = tmp_path / "sequences.jsonl"
    seeds = _load_real_seeds()
    aggregate(
        per_design_fastas=[fasta],
        design_records=[
            {
                "design_id": "design_00000",
                "backbone_pdb": "x.pdb",
                "spliced_pdb": "y.pdb",
                "binder_length": 3,
            }
        ],
        cycle=2,
        seeds=seeds,
        out_jsonl=out,
    )
    records = [json.loads(line) for line in out.read_text().splitlines()]
    assert len(records) == 1
    assert records[0]["seq"] == "EEE"
    assert records[0]["mpnn_global_score"] == pytest.approx(0.60)


def test_aggregate_emits_optional_header_fields(tmp_path: Path) -> None:
    fasta = tmp_path / "design_00000.fa"
    fasta.write_text(
        ">scaffold_input\nGGGGGG\n"
        ">T=0.1, sample=3, score=0.75, global_score=0.80, seq_recovery=0.55\nMAEELA\n"
    )
    out = tmp_path / "sequences.jsonl"
    seeds = _load_real_seeds()
    aggregate(
        per_design_fastas=[fasta],
        design_records=[
            {
                "design_id": "design_00000",
                "backbone_pdb": "x.pdb",
                "spliced_pdb": "y.pdb",
            }
        ],
        cycle=2,
        seeds=seeds,
        out_jsonl=out,
    )
    records = [json.loads(line) for line in out.read_text().splitlines()]
    assert len(records) == 1
    rec = records[0]
    assert rec["mpnn_nll"] == pytest.approx(0.75)
    assert rec["mpnn_global_score"] == pytest.approx(0.80)
    assert rec["mpnn_sample"] == 3
    assert rec["mpnn_seq_recovery"] == pytest.approx(0.55)


def test_aggregate_skips_scaffold_record(tmp_path: Path) -> None:
    fasta = tmp_path / "design_00000.fa"
    fasta.write_text(">input_scaffold\nGGGG\n" ">T=0.1, sample=1, score=0.50\nMAEELA\n")
    out = tmp_path / "sequences.jsonl"
    seeds = _load_real_seeds()
    aggregate(
        per_design_fastas=[fasta],
        design_records=[
            {
                "design_id": "design_00000",
                "backbone_pdb": "x.pdb",
                "spliced_pdb": "y.pdb",
            }
        ],
        cycle=2,
        seeds=seeds,
        out_jsonl=out,
    )
    records = [json.loads(line) for line in out.read_text().splitlines()]
    assert len(records) == 1
    assert records[0]["seq"] == "MAEELA"
