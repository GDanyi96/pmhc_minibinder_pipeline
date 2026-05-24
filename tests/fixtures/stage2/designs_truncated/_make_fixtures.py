"""Deterministic fixture generator for Stage 2 sub-run A (truncated) mock tests.

Cycle 03 sub-run A counterpart to ``../designs/_make_fixtures.py``. The
truncated in-silico funnel keeps the BAKER groove-only geometry throughout:
RFdiffusion emits 3-chain designs (A=binder, B=HLA[1:180], C=peptide), splice
re-composes them as A=HLA, B=peptide, C=binder, and ColabFold returns that same
3-chain layout (LAYOUT_CHAINS["truncated"], binder=C). beta2m is never
reintroduced in silico.

Produces, under tests/fixtures/stage2/designs_truncated/:
- mock_truncated_pmhc.pdb (cleaned BAKER target: B=HLA 4 res, C=peptide 2 res)
- mock_truncated_target.yaml (TargetManifest pointing at mock_truncated_pmhc.pdb)
- stage1_backbones/design_NNNNN.pdb (3-chain A=binder, B=HLA, C=peptide)
- subrun_summary.json (no halt_rule; fraction_geometry_pass drives the verdict
  via run_stage2._stage1_verdict, mirroring run_stage1_subrun output)
- stage1_designs.jsonl (named to NOT match run_stage2's designs.jsonl glob, so
  the mock path trusts the backbone glob -- same convention as ../designs/)
- mpnn_outputs/design_NNNNN.fa (1 scaffold record + N sampled sequences)
- af2_predictions/{design_id}_seq{K}/ (3-chain A=HLA,B=peptide,C=binder PDB +
  scores.json). All predictions pass the recalibrated truncated gate
  (iPAE<6.0 AND ipLDDT>92.0).

Run from repo root:
    uv run python tests/fixtures/stage2/designs_truncated/_make_fixtures.py
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import yaml

FIXTURE_DIR = Path(__file__).resolve().parent
N_BACKBONES = 3
N_SEQS_PER_BACKBONE = 2

# Recalibrated cycle-03 truncated gate is iPAE<6.0 AND ipLDDT>92.0
# (configs/af2_stage2.yaml). Keep every mock prediction comfortably passing.
PASS_IPAE = 3.0
PASS_IPLDDT = 95.0


def _atom_line(serial: int, chain: str, resnum: int, xyz: tuple[float, float, float]) -> str:
    x, y, z = xyz
    return (
        f"ATOM  {serial:>5d}  CA  ALA {chain}{resnum:>4d}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C"
    )


def _write_pdb(path: Path, atoms: list[tuple[str, int, tuple[float, float, float]]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["HEADER    STAGE2 SUBRUN-A TRUNCATED MOCK FIXTURE"]
    for serial, (chain, resnum, xyz) in enumerate(atoms, start=1):
        lines.append(_atom_line(serial, chain, resnum, xyz))
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def write_mock_truncated_pmhc(path: Path) -> None:
    """Cleaned BAKER groove-only target: B=HLA (4 res), C=peptide (2 res)."""
    atoms = [
        ("B", 65, (0.0, 0.0, 0.0)),
        ("B", 66, (3.8, 0.0, 0.0)),
        ("B", 150, (7.6, 0.0, 0.0)),
        ("B", 155, (11.4, 0.0, 0.0)),
        ("C", 1, (0.0, 0.0, 5.0)),
        ("C", 4, (3.0, 3.0, 5.0)),
    ]
    _write_pdb(path, atoms)


def write_mock_truncated_target(path: Path) -> None:
    """Minimal TargetManifest pointing at mock_truncated_pmhc.pdb.

    Stage 2 only reads ``primary.target_id`` and ``primary.cleaned_pdb``; the
    chains block is metadata for schema validity. The truncated splice/FASTA
    path reads chains B (HLA) and C (peptide) from the cleaned PDB itself.
    """
    primary = {
        "target_id": "wt1_a0201",
        "pdb_id": "MOCK",
        "cleaned_pdb": "mock_truncated_pmhc.pdb",
        "peptide_sequence": "RM",
        "chains": {
            "A": {"role": "heavy_chain", "length": 4},
            "B": {"role": "beta2m", "length": 2},
            "C": {"role": "peptide", "length": 2},
            "D": {"role": "binder", "length": None},
        },
        "hotspots": {
            "peptide": [{"token": "C1", "chain": "C", "resnum_clean": 1, "resnum_rcsb": 1}],
            "mhc": [{"token": "B65", "chain": "B", "resnum_clean": 65, "resnum_rcsb": 65}],
        },
        "binder": {"length_min": 1, "length_max": 10},
    }
    manifest: dict[str, Any] = {
        "cycle": "99",
        "mock": True,
        "primary": primary,
        "positive_control": dict(primary, target_id="ny1b04"),
    }
    path.write_text(yaml.safe_dump(manifest, sort_keys=True))


def write_stage1_backbones(out_dir: Path) -> list[Path]:
    """3-chain sub-run A designs: A=binder, B=HLA, C=peptide.

    The splice (splice_binder_subrun_a) takes the binder from chain A and the
    HLA/peptide from the canonical truncated target, so the design's B/C are
    intentionally small placeholders -- exercising the "design B/C ignored,
    sourced from cleaned target" contract.
    """
    hla_peptide = [
        ("B", 65, (0.0, 0.0, 0.0)),
        ("B", 66, (3.8, 0.0, 0.0)),
        ("C", 1, (0.0, 0.0, 5.0)),
    ]
    paths: list[Path] = []
    for i in range(N_BACKBONES):
        resnums = list(range(1, 4 + i % 3))  # 4..6 residues
        binder_atoms = [("A", n, (10.0 + n * 0.1, 5.0 + i * 0.2, 0.0)) for n in resnums]
        path = out_dir / f"design_{i:05d}.pdb"
        _write_pdb(path, binder_atoms + hla_peptide)
        paths.append(path)
    return paths


def write_subrun_summary(out_path: Path, designs_jsonl: Path) -> None:
    """Sub-run A summary: no halt_rule (mirrors run_stage1_subrun output).

    run_stage2._stage1_verdict derives PASS from fraction_geometry_pass >=
    Stage 1 HALT_THRESHOLD (0.10).
    """
    summary = {
        "cycle": 99,
        "subrun": "a",
        "target_id": "wt1_a0201",
        "n_attempted": N_BACKBONES,
        "n_completed": N_BACKBONES,
        "n_geometry_pass": N_BACKBONES,
        "fraction_geometry_pass": 1.0,
        "length_dist": {"mean": 5.0, "min": 4, "max": 6},
        "partial_T": 15,
        "seed_offset": 0,
        "rfdiffusion_git_sha": "mock",
        "weights_sha256": {"mock_ckpt.pt": "0" * 64},
        "wall_minutes_total": 0.01,
        "mock": True,
    }
    out_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    with designs_jsonl.open("w") as fh:
        for i in range(N_BACKBONES):
            record = {
                "design_id": f"design_{i:05d}",
                "pdb_path": f"designs/design_{i:05d}.pdb",
                "subrun": "a",
                "geometry_pass": True,
                "binder_length": 4 + i % 3,
            }
            fh.write(json.dumps(record, sort_keys=True) + "\n")


def write_mpnn_outputs(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for i in range(N_BACKBONES):
        fasta = out_dir / f"design_{i:05d}.fa"
        lines = [">scaffold_input", "GGGGGGGGGG"]
        for k in range(N_SEQS_PER_BACKBONE):
            nll = 0.50 + (i * 0.01) + (k * 0.001)
            lines.append(
                f">T=0.1, sample={k + 1}, score={nll:.4f}, "
                "global_score=0.852, model_name=v_48_020"
            )
            lines.append("MAEELRKLAE")
        fasta.write_text("\n".join(lines) + "\n")


def _build_pae_matrix(ipae_value: float, chain_lengths: dict[str, int]) -> list[list[float]]:
    """Square PAE: within-chain 0.5, binder(C)<->target(A/B) at ipae_value."""
    sizes = [chain_lengths[c] for c in ("A", "B", "C")]
    total = sum(sizes)
    matrix = [[0.5 for _ in range(total)] for _ in range(total)]
    c_start = sum(sizes[:2])
    for i in range(c_start, total):
        for j in range(0, c_start):
            matrix[i][j] = ipae_value
            matrix[j][i] = ipae_value
    return matrix


def _build_plddt(iplddt_value: float, chain_lengths: dict[str, int]) -> list[float]:
    plddt: list[float] = []
    for cid in ("A", "B"):
        plddt.extend([95.0] * chain_lengths[cid])
    plddt.extend([iplddt_value] * chain_lengths["C"])
    return plddt


def write_af2_predictions(out_dir: Path) -> None:
    """Mock prediction subdirs, 3-chain A=HLA, B=peptide, C=binder."""
    chain_lengths = {"A": 4, "B": 2, "C": 3}
    pred_atoms = [
        ("A", 1, (0.0, 0.0, 0.0)),
        ("A", 2, (3.8, 0.0, 0.0)),
        ("A", 3, (7.6, 0.0, 0.0)),
        ("A", 4, (11.4, 0.0, 0.0)),
        ("B", 1, (0.0, 0.0, 5.0)),
        ("B", 2, (3.0, 3.0, 5.0)),
        ("C", 1, (2.0, 0.0, 2.0)),
        ("C", 2, (2.0, 3.0, 2.0)),
        ("C", 3, (2.0, 6.0, 2.0)),
    ]
    for design_index in range(N_BACKBONES):
        for seq_index in range(N_SEQS_PER_BACKBONE):
            pred_id = f"design_{design_index:05d}_seq{seq_index}"
            pred_dir = out_dir / pred_id
            pred_dir.mkdir(parents=True, exist_ok=True)
            scores = {
                "pae": _build_pae_matrix(PASS_IPAE, chain_lengths),
                "plddt": _build_plddt(PASS_IPLDDT, chain_lengths),
            }
            (
                pred_dir / f"{pred_id}_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json"
            ).write_text(json.dumps(scores))
            _write_pdb(
                pred_dir
                / f"{pred_id}_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb",
                pred_atoms,
            )


def main() -> None:
    write_mock_truncated_pmhc(FIXTURE_DIR / "mock_truncated_pmhc.pdb")
    write_mock_truncated_target(FIXTURE_DIR / "mock_truncated_target.yaml")
    write_stage1_backbones(FIXTURE_DIR / "stage1_backbones")
    write_subrun_summary(
        FIXTURE_DIR / "subrun_summary.json",
        FIXTURE_DIR / "stage1_designs.jsonl",
    )
    write_mpnn_outputs(FIXTURE_DIR / "mpnn_outputs")
    write_af2_predictions(FIXTURE_DIR / "af2_predictions")
    print(f"wrote stage2 sub-run A truncated fixtures under {FIXTURE_DIR}")


if __name__ == "__main__":
    main()
