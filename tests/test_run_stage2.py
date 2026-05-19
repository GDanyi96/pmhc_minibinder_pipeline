"""Unit tests for the new real-mode helpers in scripts/run_stage2.py.

These cover the pure-Python pieces (FASTA writer, chain-sequence
extraction, binder-length manifest loader) without touching the
ProteinMPNN / ColabFold subprocesses, which only exist on the pod.
"""

from __future__ import annotations

import json
from pathlib import Path

from scripts.run_stage2 import (
    _build_multimer_fasta,
    _extract_target_chain_sequences,
    _load_binder_lengths,
)


def test_extract_target_chain_sequences_from_mock_pmhc() -> None:
    """The mock pMHC fixture has 2 ALA residues per chain across A/B/C."""
    pdb = Path("tests/fixtures/stage2/designs/mock_pmhc.pdb")
    seqs = _extract_target_chain_sequences(pdb)
    assert seqs == {"A": "AA", "B": "AA", "C": "AA"}


def test_build_multimer_fasta_colon_separated_single_record(tmp_path: Path) -> None:
    """Trap #7: ColabFold multimer FASTA is colon-separated within a
    single record. Chain order A:B:C:D, no trailing colon, one >header
    per file.
    """
    target_seqs = {"A": "AAAA", "B": "BB", "C": "C", "D": "ignored"}
    fasta = tmp_path / "design_00042_seq01.fasta"
    _build_multimer_fasta(fasta, "design_00042_seq01", target_seqs, "DDDDD")
    content = fasta.read_text()
    assert content == ">design_00042_seq01\nAAAA:BB:C:DDDDD\n"
    assert content.count(">") == 1
    assert not content.rstrip().endswith(":")


def test_build_multimer_fasta_uses_pdb_chain_order(tmp_path: Path) -> None:
    """Chain order must follow PDB convention A->B->C->D regardless of
    dict insertion order. If the helper trusts dict insertion, a target
    dict ordered C/B/A would silently break the AF2 input.
    """
    target_seqs = {"C": "CCC", "B": "BB", "A": "A"}
    fasta = tmp_path / "rec.fasta"
    _build_multimer_fasta(fasta, "rec", target_seqs, "DDDD")
    body = fasta.read_text().splitlines()[1]
    chains = body.split(":")
    assert chains == ["A", "BB", "CCC", "DDDD"]


def test_load_binder_lengths_from_designs_jsonl(tmp_path: Path) -> None:
    jsonl = tmp_path / "designs.jsonl"
    jsonl.write_text(
        json.dumps({"design_id": "design_00000", "binder_length": 78})
        + "\n"
        + json.dumps({"design_id": "design_00001", "binder_length": 85})
        + "\n"
    )
    lengths = _load_binder_lengths(jsonl)
    assert lengths == {"design_00000": 78, "design_00001": 85}


def test_load_binder_lengths_ignores_blank_lines(tmp_path: Path) -> None:
    jsonl = tmp_path / "designs.jsonl"
    jsonl.write_text("\n" + json.dumps({"design_id": "design_00000", "binder_length": 78}) + "\n\n")
    assert _load_binder_lengths(jsonl) == {"design_00000": 78}
