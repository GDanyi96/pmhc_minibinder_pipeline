"""Tests for the cycle-03 Stage 2 sub-run A truncated path (3-chain designs).

Mirrors the cycle-02 Stage 2 mock canary but for the BAKER groove-only layout:
RFdiffusion 3-chain output (A=binder, B=HLA, C=peptide) -> splice to A=HLA,
B=peptide, C=binder -> ProteinMPNN designs chain C -> AF2 FASTA HLA:peptide:binder
-> compute_metrics with target_layout="truncated" (binder=C, peptide=B, mhc=A).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from scripts import run_stage2
from workflow.scripts import compute_metrics, splice_binder

REPO_ROOT = Path(__file__).resolve().parent.parent
TRUNC_FIXTURES = REPO_ROOT / "tests/fixtures/stage2/designs_truncated"


def _write_ca_pdb(path: Path, chain_lengths: dict[str, int]) -> None:
    lines: list[str] = []
    serial = 1
    for chain_id, n in chain_lengths.items():
        for resnum in range(1, n + 1):
            x, y, z = float(serial), 0.0, 0.0
            lines.append(
                f"ATOM  {serial:>5d}  CA  ALA {chain_id}{resnum:>4d}    "
                f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C"
            )
            serial += 1
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


# --- _stage1_verdict adapter --------------------------------------------------


def test_stage1_verdict_uses_explicit_halt_rule() -> None:
    assert run_stage2._stage1_verdict({"halt_rule": {"verdict": "PASS"}}) == "PASS"
    assert run_stage2._stage1_verdict({"halt_rule": {"verdict": "FAIL"}}) == "FAIL"


def test_stage1_verdict_derives_from_fraction_for_subrun_summary() -> None:
    # Sub-run summaries carry fraction_geometry_pass but no halt_rule.
    assert run_stage2.STAGE1_HALT_THRESHOLD == 0.10
    assert run_stage2._stage1_verdict({"fraction_geometry_pass": 0.10}) == "PASS"
    assert run_stage2._stage1_verdict({"fraction_geometry_pass": 0.48}) == "PASS"
    assert run_stage2._stage1_verdict({"fraction_geometry_pass": 0.05}) == "FAIL"


# --- splice_binder_subrun_a ---------------------------------------------------


def test_splice_subrun_a_reorders_to_hla_peptide_binder(tmp_path: Path) -> None:
    """3-chain design A=binder(5), B=HLA(7), C=peptide(3) + canonical truncated
    target B=HLA(7), C=peptide(3) -> spliced A=HLA, B=peptide, C=binder."""
    design = tmp_path / "design_3000.pdb"
    _write_ca_pdb(design, {"A": 5, "B": 7, "C": 3})  # design B/C are placeholders
    target = tmp_path / "trunc_target.pdb"
    _write_ca_pdb(target, {"B": 7, "C": 3})  # canonical HLA/peptide
    out = tmp_path / "spliced.pdb"

    splice_binder.splice_binder_subrun_a(design, target, out)

    from Bio.PDB import PDBParser

    model = PDBParser(QUIET=True).get_structure("s", str(out))[0]
    lengths = {c.id: len([r for r in c.get_residues() if r.id[0].strip() == ""]) for c in model}
    assert lengths == {"A": 7, "B": 3, "C": 5}  # HLA, peptide, binder


def test_splice_subrun_a_requires_binder_on_chain_a(tmp_path: Path) -> None:
    design = tmp_path / "no_a.pdb"
    _write_ca_pdb(design, {"B": 7, "C": 3})
    target = tmp_path / "trunc_target.pdb"
    _write_ca_pdb(target, {"B": 7, "C": 3})
    with pytest.raises(ValueError, match="binder on chain 'A'"):
        splice_binder.splice_binder_subrun_a(design, target, tmp_path / "out.pdb")


def test_splice_subrun_a_rejects_four_chain_design(tmp_path: Path) -> None:
    """A cycle-02 4-chain design must not silently flow through the truncated splice."""
    design = tmp_path / "four_chain.pdb"
    _write_ca_pdb(design, {"A": 5, "B": 7, "C": 3, "D": 9})
    target = tmp_path / "trunc_target.pdb"
    _write_ca_pdb(target, {"B": 7, "C": 3})
    with pytest.raises(ValueError, match="unexpected chain"):
        splice_binder.splice_binder_subrun_a(design, target, tmp_path / "out.pdb")


# --- FASTA chain order --------------------------------------------------------


def test_truncated_fasta_is_hla_peptide_binder(tmp_path: Path) -> None:
    target_seqs = {"B": "HHHHHH", "C": "PP"}  # HLA, peptide
    fasta = tmp_path / "rec.fasta"
    run_stage2._build_multimer_fasta(fasta, "rec", target_seqs, "BBBBB", target_layout="truncated")
    body = fasta.read_text().splitlines()[1]
    chains = body.split(":")
    assert chains == ["HHHHHH", "PP", "BBBBB"]  # HLA : peptide : binder
    assert len(chains) == 3


def test_full_fasta_unchanged_default() -> None:
    # The default (full) layout still produces the cycle-02 A:B:C:D order.
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        fasta = Path(d) / "rec.fasta"
        run_stage2._build_multimer_fasta(fasta, "rec", {"A": "AA", "B": "BB", "C": "C"}, "DDDD")
        assert fasta.read_text().splitlines()[1] == "AA:BB:C:DDDD"


# --- MPNN binder-chain dispatch ----------------------------------------------


def test_mpnn_assign_uses_chain_c_for_truncated(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    spliced = tmp_path / "design_00000.pdb"
    _write_ca_pdb(spliced, {"A": 7, "B": 3, "C": 5})
    captured: list[list[str]] = []

    monkeypatch.setattr(run_stage2, "_resolve_proteinmpnn_python", lambda: "python")
    monkeypatch.setattr(run_stage2, "_stage_input_pdbs", lambda pdbs, wd: tmp_path / "input")

    def fake_run(cmd, **kwargs):  # type: ignore[no-untyped-def]
        captured.append(cmd)

    monkeypatch.setattr(run_stage2.subprocess, "run", fake_run)
    # Stub the FASTA collection: write the expected per-design .fa so the
    # collector at the end of _run_real_proteinmpnn finds it.
    seqs_dir = tmp_path / "mpnn_work" / "mpnn_out" / "seqs"
    seqs_dir.mkdir(parents=True)
    (seqs_dir / "design_00000.fa").write_text(">x\nAAAA\n")

    run_stage2._run_real_proteinmpnn(
        [spliced],
        tmp_path / "mpnn_work",
        {"n_seqs_per_backbone": 2, "sampling_temperature": 0.1},
        seed=3200,
        binder_chain="C",
    )
    assign = next(c for c in captured if any("assign_fixed_chains.py" in a for a in c))
    idx = assign.index("--chain_list")
    assert assign[idx + 1] == "C"


# --- metrics decomposition (real path, truncated layout) ----------------------


def test_metrics_for_design_truncated_decomposition() -> None:
    pred_dir = TRUNC_FIXTURES / "af2_predictions" / "design_00000_seq0"
    m = compute_metrics.metrics_for_design(pred_dir, target_layout="truncated")
    for key in ("ipae", "iplddt", "ppi_pae_int_peptide", "ppi_pae_int_mhc"):
        assert key in m
        assert not math.isnan(float(m[key]))
        assert math.isfinite(float(m[key]))
    # Binder (C) is in contact with both peptide (B) and HLA (A) in the fixture.
    assert m["ipae"] < 6.0


# --- mock end-to-end ----------------------------------------------------------


def test_subrun_a_stage2_mock_end_to_end(tmp_path: Path) -> None:
    rc = run_stage2.run(
        stage1_summary_path=TRUNC_FIXTURES / "subrun_summary.json",
        target_manifest=TRUNC_FIXTURES / "mock_truncated_target.yaml",
        mpnn_config=REPO_ROOT / "configs/proteinmpnn_cycle03.yaml",
        af2_config=REPO_ROOT / "configs/af2_stage2.yaml",
        seeds_config=REPO_ROOT / "configs/seeds.yaml",
        cycle=99,
        out_dir=tmp_path / "stage2_subrun_a",
        mock=True,
        metrics_only=False,
        target_layout="truncated",
    )
    assert rc == 0

    # Spliced PDBs must be 3-chain A=HLA, B=peptide, C=binder.
    from Bio.PDB import PDBParser

    spliced = sorted(
        (tmp_path / "stage2_subrun_a" / "proteinmpnn_designs" / "spliced").glob("design_*.pdb")
    )
    assert spliced
    model = PDBParser(QUIET=True).get_structure("s", str(spliced[0]))[0]
    assert {c.id for c in model.get_chains()} == {"A", "B", "C"}

    # Metrics must be non-NaN and computed under the truncated layout (binder=C).
    # If "full" had leaked through, binder=D would be absent in a 3-chain PDB and
    # iPAE would be +inf.
    metrics = [
        json.loads(line)
        for line in (tmp_path / "stage2_subrun_a" / "af2_designs" / "metrics.jsonl")
        .read_text()
        .splitlines()
        if line.strip()
    ]
    assert metrics
    for rec in metrics:
        assert math.isfinite(float(rec["iPAE"]))
        assert math.isfinite(float(rec["ipLDDT"]))
        assert float(rec["iPAE"]) < 6.0

    summary = json.loads((tmp_path / "stage2_subrun_a" / "stage2_summary.json").read_text())
    assert summary["halt_rule"]["verdict"] == "PASS"
