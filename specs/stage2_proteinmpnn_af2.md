# Stage 2 — Sequence Design (ProteinMPNN) + Structure Validation (AF2-multimer)

## Goal

Design sequences on each RFdiffusion backbone using ProteinMPNN in complex
mode (fixing MHC chains), then validate the binder–pMHC complex with
AF2-multimer via ColabFold. Filter by iPAE, ipLDDT, and ProteinMPNN NLL.
Include the P1–N2 controls panel on every cycle.

This stage has three sub-rules: `02_proteinmpnn`, `03_colabfold`, `04_metrics`.

## Inputs

- `results/cycle_NN/stage1/backbones/design_{i:05d}.pdb` — RFdiffusion outputs.
- `results/cycle_NN/stage0/target.yaml` — chain IDs, hotspots.
- `data/controls/controls_manifest.yaml` — P1–N2 binder sequences.
- `configs/thresholds.yaml` — iPAE, ipLDDT, NLL cuts.

## Outputs

- `results/cycle_NN/stage2/proteinmpnn/design_{i:05d}.fasta` — 8 sequences/backbone.
- `results/cycle_NN/stage2/colabfold/design_{i:05d}_seq{j}/ranked_0.pdb`
- `results/cycle_NN/stage2/colabfold/design_{i:05d}_seq{j}/scores.json` —
  pae, plddt, ranking confidence.
- `results/cycle_NN/stage2/metrics.parquet` — per-design metrics table.
- `results/cycle_NN/stage2/survivors.json` — designs passing all filters.
- `results/cycle_NN/stage2/controls_metrics.json` — P1–N2 outcomes.

## Tools

- ProteinMPNN: upstream `dauparas/ProteinMPNN` cloned natively on the pod
  (no DinD); invoked via `python $PROTEINMPNN_DIR/protein_mpnn_run.py …`.
  Complex mode: `--fix_chains "A B C"` (only chain D = binder is designed).
- ColabFold installed natively via the `colabfold` optional extra
  (`colabfold_batch` CLI); remote MMseqs2 MSA. AF2-multimer preset,
  `num_recycles=6`.
- `pandas` / `pyarrow` for parquet output.

## Anchored references

- `HADRUP_JENKINS_2025` — final cuts iPAE < 7 Å, ipLDDT > 92; T=0.1 for MPNN.
- `MARES_IOANNIDIS_2025` — AF2 `num_recycles=6`, multimer preset.
- `BAKER_LAB_2025` — sequences per backbone = 8 (then top-K by NLL).

## Filter thresholds (cycle 1)

| Metric           | Cycle 1 cut | Final cut       |
|------------------|-------------|-----------------|
| AF2 iPAE         | < 12 Å      | **< 7 Å**       |
| AF2 ipLDDT       | loose       | **> 92**        |
| MPNN NLL         | top 10 %    | top 5 %         |

All values from `configs/thresholds.yaml`. Do not hard-code.

## Controls panel (mandatory)

| #  | Binder                                  | Target           | Expected                |
|----|-----------------------------------------|------------------|-------------------------|
| P1 | Baker WT1 binder                        | WT1/A\*02:01     | iPAE<7, top 10 %         |
| P2 | Jenkins NY1-B04                         | SLLMWITQC/A\*02:01 | iPAE ≈ 6.5             |
| P3 | Baker WT1 binder                        | MART-1/A\*02:01  | iPAE > P1 + 3 Å         |
| N1 | Scrambled Baker WT1                     | WT1/A\*02:01     | iPAE > 15               |
| N2 | Random 65-aa (natural AA freqs)         | WT1/A\*02:01     | iPAE > 20               |

If P1's iPAE rank falls below 50th percentile of the design population, the
pipeline halts and `results/cycle_NN/HALT.txt` is written. Stage 6 reports the halt.

## Implementation tasks

- [ ] `run_proteinmpnn.py`: glue script around ProteinMPNN CLI; pass
      `--fixed_chains "A B C"`, `T=0.1`, generate 8 sequences/backbone.
- [ ] `run_colabfold.py`: invoke ColabFold per (backbone, sequence); use
      remote MSA (`--msa-mode mmseqs2_uniref_env`); `--num-recycle 6`.
- [ ] `compute_metrics.py`: parse ColabFold scores.json, compute iPAE
      (interface PAE between chain D and chains A/B/C), ipLDDT (pLDDT
      averaged over interface residues), combine with MPNN NLL.
- [ ] Apply cycle-1 filters; write survivors.json and metrics.parquet.
- [ ] Run controls in parallel; emit controls_metrics.json.
- [ ] Implement halt rule (P1 percentile check).
- [ ] `--mock` flags on all three scripts: cp fixtures to outputs.

## Verification criteria

- Smoke test (each script with `--mock`): exits 0 in <1s, produces stub
  outputs matching expected paths.
- Real test (pod, single backbone): one full design → AF2 cycle completes
  in <10 min; metrics row populated.
- Controls test: cycle 1 on real designs has P1 iPAE < 7 (sanity check).

## Pitfalls

- AF2 chain order: `A=HC, B=β2m, C=peptide, D=binder` — wrong order silently
  yields garbage iPAE.
- ColabFold remote MSA can rate-limit at >50 concurrent jobs. Throttle.
- ProteinMPNN's "complex mode" needs the target re-attached to the backbone
  PDB (RFdiffusion outputs binder only by default).
- iPAE definition varies across papers. Use `HADRUP_JENKINS_2025` Fig 1B:
  mean PAE over chain-D residues to chain-{A,B,C} residues, both ways.
- ipLDDT must be averaged over **interface** residues only (binder residues
  within 8 Å of MHC heavy atoms), not the entire binder.
