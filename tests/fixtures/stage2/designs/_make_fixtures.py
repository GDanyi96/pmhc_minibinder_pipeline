"""Deterministic fixture generator for Stage 2 real-designs mock tests.

Produces, under tests/fixtures/stage2/designs/:
- stage1_backbones/design_NNNNN.pdb (10 single-chain binder PDBs, chain "A")
- stage1_summary.json (verdict=PASS, n_completed=10, fraction_geometry_pass=1.0)
- stage1_designs.jsonl (10 minimal records)
- mock_pmhc.pdb (cleaned pMHC: A=2 residues HC, B=2 beta2m, C=2 peptide)
- mock_target.yaml (TargetManifest pointing at mock_pmhc.pdb)
- mpnn_outputs/design_NNNNN.fa (10 FASTAs, each with 1 scaffold + 4 sampled
  sequences, deterministic NLLs)
- af2_predictions/{design_id}_seq{K}/ (40 ColabFold-shaped subdirs;
  each contains a synthetic scores.json + a minimal 4-chain PDB)

Calibration: exactly 4 of 40 mock predictions pass iPAE<12 AND ipLDDT>88
-> observed = 0.10 (boundary). See the comment block in
build_mock_metrics() below.

Run from repo root:
    uv run python tests/fixtures/stage2/designs/_make_fixtures.py
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import yaml

FIXTURE_DIR = Path(__file__).resolve().parent
N_BACKBONES = 10
N_SEQS_PER_BACKBONE = 4
N_TOTAL_RECORDS = N_BACKBONES * N_SEQS_PER_BACKBONE  # 40

# Boundary-test fixture: exactly 4 of 40 mock predictions pass the intermediate
# cut (iPAE < halt_cut_ipae_max AND ipLDDT > halt_cut_iplddt_min from
# configs/af2_stage2.yaml) -> observed = 0.10 exactly. This tests the halt rule
# uses >= (PASS at boundary), not > (FAIL at boundary). Any refactor that
# computes `observed` differently (e.g. integer arithmetic, percentage scaling)
# must preserve `observed == 0.10` and the corresponding PASS verdict. Do not
# change pass count without updating the halt-rule semantic test in
# tests/test_stage2_designs_halt_gate.py. PASS_* are kept comfortably inside the
# cycle-03 recalibrated gate (6.0 / 92.0); FAIL_* comfortably outside.
N_PASSING = 4
PASS_IPAE = 3.0
PASS_IPLDDT = 95.0
FAIL_IPAE = 18.0
FAIL_IPLDDT = 60.0


def _atom_line(serial: int, chain: str, resnum: int, xyz: tuple[float, float, float]) -> str:
    x, y, z = xyz
    return (
        f"ATOM  {serial:>5d}  CA  ALA {chain}{resnum:>4d}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C"
    )


def _write_pdb(path: Path, atoms: list[tuple[str, int, tuple[float, float, float]]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["HEADER    STAGE2 DESIGNS MOCK FIXTURE"]
    for serial, (chain, resnum, xyz) in enumerate(atoms, start=1):
        lines.append(_atom_line(serial, chain, resnum, xyz))
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def write_mock_pmhc(path: Path) -> None:
    """Cleaned pMHC fixture: A=2 (HC), B=2 (beta2m), C=2 (peptide)."""
    atoms = [
        ("A", 65, (0.0, 0.0, 0.0)),
        ("A", 66, (3.8, 0.0, 0.0)),
        ("B", 1, (20.0, 0.0, 0.0)),
        ("B", 2, (20.0, 3.8, 0.0)),
        ("C", 1, (0.0, 0.0, 5.0)),
        ("C", 4, (3.0, 3.0, 5.0)),
    ]
    _write_pdb(path, atoms)


def write_mock_target(path: Path) -> None:
    """Minimal TargetManifest pointing at mock_pmhc.pdb (relative path)."""
    manifest: dict[str, Any] = {
        "cycle": "99",
        "mock": True,
        "primary": {
            "target_id": "wt1_a0201",
            "pdb_id": "MOCK",
            "cleaned_pdb": "mock_pmhc.pdb",
            "peptide_sequence": "RM",
            "chains": {
                "A": {"role": "heavy_chain", "length": 2},
                "B": {"role": "beta2m", "length": 2},
                "C": {"role": "peptide", "length": 2},
                "D": {"role": "binder", "length": None},
            },
            "hotspots": {
                "peptide": [{"token": "C1", "chain": "C", "resnum_clean": 1, "resnum_rcsb": 1}],
                "mhc": [{"token": "A65", "chain": "A", "resnum_clean": 65, "resnum_rcsb": 65}],
            },
            "binder": {"length_min": 1, "length_max": 10},
        },
        "positive_control": {
            "target_id": "ny1b04",
            "pdb_id": "MOCK",
            "cleaned_pdb": "mock_pmhc.pdb",
            "peptide_sequence": "RM",
            "chains": {
                "A": {"role": "heavy_chain", "length": 2},
                "B": {"role": "beta2m", "length": 2},
                "C": {"role": "peptide", "length": 2},
                "D": {"role": "binder", "length": None},
            },
            "hotspots": {
                "peptide": [{"token": "C1", "chain": "C", "resnum_clean": 1, "resnum_rcsb": 1}],
                "mhc": [{"token": "A65", "chain": "A", "resnum_clean": 65, "resnum_rcsb": 65}],
            },
            "binder": {"length_min": 1, "length_max": 10},
        },
    }
    path.write_text(yaml.safe_dump(manifest, sort_keys=True))


def write_stage1_backbones(out_dir: Path) -> list[Path]:
    """10 four-chain (A=HC, B=beta2m, C=peptide, D=binder) PDBs.

    Mirrors the cycle 02 RFdiffusion output layout: the motif chains
    A/B/C are passed through unchanged, and the designed binder is
    written to chain D. A/B/C atom counts match the cleaned pMHC
    fixture so the splicer's "A/B/C ignored, sourced from cleaned"
    contract is exercised end-to-end (Trap #29).
    """
    pmhc_atoms: list[tuple[str, int, tuple[float, float, float]]] = [
        ("A", 65, (0.0, 0.0, 0.0)),
        ("A", 66, (3.8, 0.0, 0.0)),
        ("B", 1, (20.0, 0.0, 0.0)),
        ("B", 2, (20.0, 3.8, 0.0)),
        ("C", 1, (0.0, 0.0, 5.0)),
        ("C", 4, (3.0, 3.0, 5.0)),
    ]
    paths: list[Path] = []
    for i in range(N_BACKBONES):
        resnums = list(range(1, 3 + i % 3))  # 3..5 residues
        binder_atoms: list[tuple[str, int, tuple[float, float, float]]] = [
            ("D", n, (10.0 + n * 0.1, 5.0 + i * 0.2, 0.0)) for n in resnums
        ]
        path = out_dir / f"design_{i:05d}.pdb"
        _write_pdb(path, pmhc_atoms + binder_atoms)
        paths.append(path)
    return paths


def write_stage1_summary(out_path: Path, designs_jsonl: Path) -> None:
    summary = {
        "cycle": 99,
        "target_id": "wt1_a0201",
        "n_attempted": N_BACKBONES,
        "n_completed": N_BACKBONES,
        "n_geometry_pass": N_BACKBONES,
        "fraction_geometry_pass": 1.0,
        "length_dist": {"mean": 4.0, "min": 3, "max": 5},
        "seed_range": [99000, 99099],
        "halt_rule": {
            "name": "geometry_pass_rate",
            "threshold": 0.50,
            "observed": 1.0,
            "verdict": "PASS",
            "note": "Stage 2 mock upstream fixture",
        },
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
                "geometry_pass": True,
                "binder_length": 3 + i % 3,
            }
            fh.write(json.dumps(record, sort_keys=True) + "\n")


def write_mpnn_outputs(out_dir: Path) -> None:
    """One FASTA per backbone, format mirrors ProteinMPNN real output:
    first record is the fixed scaffold, then N_SEQS_PER_BACKBONE samples
    with `T=0.1, sample=K, score=N`. Deterministic NLLs so the top-N
    ranking step downstream is reproducible.
    """
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
    """Build a square PAE matrix with within-chain entries at 0.5 and
    cross-chain entries between chain D and chains A/B/C set to ipae_value.
    All other cross-chain entries (A<->B, B<->C, etc.) are also 0.5 -- only
    the D<->ABC slice drives the iPAE computation.
    """
    sizes = [chain_lengths[c] for c in ("A", "B", "C", "D")]
    total = sum(sizes)
    matrix = [[0.5 for _ in range(total)] for _ in range(total)]
    d_start = sum(sizes[:3])
    d_stop = total
    for i in range(d_start, d_stop):
        for j in range(0, d_start):
            matrix[i][j] = ipae_value
            matrix[j][i] = ipae_value
    return matrix


def _build_plddt(iplddt_value: float, chain_lengths: dict[str, int]) -> list[float]:
    """plddt array: chain D residues at iplddt_value (drives ipLDDT),
    chains A/B/C at a high baseline (unused but realistic)."""
    plddt: list[float] = []
    for cid in ("A", "B", "C"):
        plddt.extend([95.0] * chain_lengths[cid])
    plddt.extend([iplddt_value] * chain_lengths["D"])
    return plddt


def write_af2_predictions(out_dir: Path) -> None:
    """40 mock prediction subdirs, exactly N_PASSING of which pass the
    intermediate cut (iPAE<12 AND ipLDDT>88). Each subdir contains:
    - {pred_id}_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb
    - {pred_id}_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json

    Pass/fail allocation: the first N_PASSING records by (design_index,
    seq_index) order get the passing metrics; the rest fail.
    """
    chain_lengths = {"A": 2, "B": 2, "C": 2, "D": 3}
    binder_atoms_d = [("D", n, (10.0 + n * 0.1, 5.0, 0.0)) for n in (1, 2, 3)]
    pmhc_atoms = [
        ("A", 1, (0.0, 0.0, 0.0)),
        ("A", 2, (3.8, 0.0, 0.0)),
        ("B", 1, (20.0, 0.0, 0.0)),
        ("B", 2, (20.0, 3.8, 0.0)),
        ("C", 1, (0.0, 0.0, 5.0)),
        ("C", 2, (3.0, 3.0, 5.0)),
    ]
    pred_atoms = pmhc_atoms + binder_atoms_d

    rank_idx = 0
    for design_index in range(N_BACKBONES):
        for seq_index in range(N_SEQS_PER_BACKBONE):
            pred_id = f"design_{design_index:05d}_seq{seq_index}"
            pred_dir = out_dir / pred_id
            pred_dir.mkdir(parents=True, exist_ok=True)
            passes = rank_idx < N_PASSING
            ipae = PASS_IPAE if passes else FAIL_IPAE
            iplddt = PASS_IPLDDT if passes else FAIL_IPLDDT

            scores = {
                "pae": _build_pae_matrix(ipae, chain_lengths),
                "plddt": _build_plddt(iplddt, chain_lengths),
            }
            scores_path = (
                pred_dir / f"{pred_id}_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json"
            )
            scores_path.write_text(json.dumps(scores))

            pdb_path = (
                pred_dir
                / f"{pred_id}_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb"
            )
            _write_pdb(pdb_path, pred_atoms)
            rank_idx += 1


def main() -> None:
    write_mock_pmhc(FIXTURE_DIR / "mock_pmhc.pdb")
    write_mock_target(FIXTURE_DIR / "mock_target.yaml")
    write_stage1_backbones(FIXTURE_DIR / "stage1_backbones")
    write_stage1_summary(
        FIXTURE_DIR / "stage1_summary.json",
        FIXTURE_DIR / "stage1_designs.jsonl",
    )
    write_mpnn_outputs(FIXTURE_DIR / "mpnn_outputs")
    write_af2_predictions(FIXTURE_DIR / "af2_predictions")
    print(f"wrote stage2 designs fixtures under {FIXTURE_DIR}")


if __name__ == "__main__":
    main()
