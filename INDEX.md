# INDEX.md — File map

One-line description per tracked file. New file? Add a line here.

## Top-level

- `CLAUDE.md` — Project rulebook (locked decisions, workflow rules, thresholds).
- `CLAUDE.local.md` — Personal runtime context (gitignored).
- `INDEX.md` — This file.
- `README.md` — Portfolio-facing project overview.
- `LICENSE` — MIT.
- `bootstrap.sh` — One-command RunPod setup (deps, model weights, ProteinMPNN clone, smoke test).
- `PIPELINE_STATUS.md` — Live per-stage status, CI/real-run state, locked architectural decisions.
- `pyproject.toml` — `uv`-managed Python project + pinned minimal deps.
- `Snakefile` — Top-level pipeline DAG; includes all `workflow/rules/*.smk`.
- `.gitignore` — Python + project-specific exclusions (CLAUDE.local.md, settings.local.json, results/).

## .claude/

- `.claude/settings.json` — Shared bash command allowlist (committed).
- `.claude/settings.local.json` — Personal allowlist (gitignored).
- `.claude/hooks/pre-push.sh` — Pre-push verification gate (snakemake dry-run).

## .github/

- `.github/workflows/ci.yml` — Lint (ruff, black, mypy --strict) + snakemake dry-run on push/PR.

## specs/

- `specs/stage0_target_specification.md` — Target + off-target panel definition, PDB cleaning, hotspot selection.
- `specs/stage1_rfdiffusion.md` — Backbone generation via RFdiffusion (native install, partial diffusion).
- `specs/stage2_proteinmpnn_af2.md` — Stage 2 real-designs contract: MPNN + AF2 fan-in + halt gate on Stage 1 outputs.
- `specs/stage2_controls.md` — Stage 2 controls contract (P1-N2 literature binders, cycle 1).
- `specs/stage3_crosspan.md` — In silico cross-reactivity panning against off-target peptide grid.
- `specs/stage4_diversity.md` — ESM-2 embeddings + farthest-point sampling diversity selection.
- `specs/stage5_active_learning.md` — LightGBM + GP surrogate, UCB acquisition for next cycle.
- `specs/stage6_reporting.md` — Cycle aggregation, publication figures, decision log.

## workflow/rules/

- `workflow/rules/00_target_prep.smk` — Download/clean target PDBs, emit target.yaml.
- `workflow/rules/01_rfdiffusion.smk` — Run RFdiffusion to generate N backbones.
- `workflow/rules/02_proteinmpnn.smk` — Design sequences on each backbone.
- `workflow/rules/03_colabfold.smk` — AF2-multimer prediction via ColabFold.
- `workflow/rules/04_metrics.smk` — Compute iPAE, ipLDDT, NLL filters.
- `workflow/rules/02b_proteinmpnn_designs.smk` — Stage 2 designs end-to-end pipeline (splice + MPNN + AF2 fan-in + halt).
- `workflow/rules/03b_af2_designs.smk` — Marker-only DAG anchor downstream of stage2_designs.
- `workflow/rules/05_crosspan.smk` — Cross-reactivity panel against off-targets.
- `workflow/rules/06_embedding.smk` — ESM-2 embeddings + FPS diversity selection.
- `workflow/rules/07_active_learning.smk` — Surrogate training + UCB acquisition.
- `workflow/rules/08_reporting.smk` — Aggregate to cycle report.

## workflow/scripts/

- `workflow/scripts/__init__.py` — Marks workflow.scripts as a package (mypy + imports).
- `workflow/scripts/prep_target.py` — Stage 0 entry point.
- `workflow/scripts/run_rfdiffusion.py` — Stage 1 entry point.
- `workflow/scripts/run_proteinmpnn.py` — Stage 2a: ProteinMPNN complex-mode wrapper. **DEPRECATED for real-mode** (cycle 02+); mock paths retained for cycle-1 controls. Real-mode bypasses this and invokes upstream ProteinMPNN directly from `scripts/run_stage2.py`.
- `workflow/scripts/run_colabfold.py` — Stage 2b: ColabFold (AF2-multimer) native wrapper.
- `workflow/scripts/compute_metrics.py` — Stage 2c: iPAE / ipLDDT / BSA computation.
- `workflow/scripts/splice_binder.py` — Stage 2 designs: splice Stage 1 binder (chain A) onto cleaned pMHC as chain D.
- `workflow/scripts/aggregate_mpnn_outputs.py` — Stage 2 designs: parse per-design MPNN FASTAs into a unified ranked sequences.jsonl; MPNN/AF2 seed helpers.
- `workflow/scripts/crosspan.py` — Stage 3 entry point.
- `workflow/scripts/embed_designs.py` — Stage 4 entry point.
- `workflow/scripts/active_learning.py` — Stage 5 entry point.
- `workflow/scripts/render_report.py` — Stage 6 entry point.

## scripts/

- `scripts/__init__.py` — Marks scripts as a package.
- `scripts/generate_negatives.py` — Deterministic N1/N2 binder sequences (seed=42).
- `scripts/run_controls.py` — Stage 2 Part A controls orchestrator + halt gate.
- `scripts/run_stage1.py` — Stage 1 RFdiffusion orchestrator + geometry-pass halt gate.
- `scripts/run_stage2.py` — Stage 2 designs orchestrator (splice + MPNN + AF2 fan-in + halt gate at fraction_pass_intermediate >= 0.10).

## configs/

- `configs/target_wt1_a0201.yaml` — Locked primary target spec (chains, hotspots, length range 70-110).
- `configs/target_2bnr_a0201.yaml` — Positive-control target spec (Jenkins NY1-B04 / 2BNR).
- `configs/thresholds.yaml` — All numerical filter thresholds (iPAE, ipLDDT, NLL, etc).
- `configs/rfdiffusion_default.yaml` — RFdiffusion noise scaling, num_designs, partial_T (contigmap built at runtime).
- `configs/seeds.yaml` — Stage 1 seed formula (`cycle * 1000 + design_index`) + per-cycle reserved ranges.
- `configs/proteinmpnn_chains.json` — Per-PDB chain assignment for ProteinMPNN complex mode.
- `configs/proteinmpnn_default.yaml` — Stage 2 designs MPNN sampling defaults (T=0.1, 4 seqs/backbone, chain D designed).
- `configs/af2_stage2.yaml` — Stage 2 designs AF2 fan-in config (top-100, intermediate cut iPAE<12/ipLDDT>88, halt 0.10).

## data/

- `data/targets/.gitkeep` — PDB files land here (downloaded by bootstrap.sh).
- `data/controls/controls_manifest.yaml` — P1, P2, P3, N1, N2 control definitions.

## docs/

- `docs/pod_quickstart.md` — Fresh-pod three-line sequence, debug flags, recovery from failed metrics steps.
- `docs/known_traps.md` — Empirical gotchas encountered bringing the pipeline up on RunPod (Docker, uv, ColabFold output layout, PAE key, JAX drift, LD_LIBRARY_PATH, AF2 specificity blind spot).

## tests/

- `tests/__init__.py` — Marks tests as a package.
- `tests/test_dry_run.py` — Asserts `snakemake --dry-run --config mock=true -j1` exits 0.
- `tests/test_halt_gate.py` — Stage 2 controls halt-gate unit tests.
- `tests/test_compute_metrics.py` — Stage 2 metrics computation tests.
- `tests/test_run_rfdiffusion.py` — Stage 1 helper unit tests (contigmap, seeds, geometry pass).
- `tests/test_stage1_halt_gate.py` — Stage 1 orchestrator + halt-gate tests (pass + fail paths).
- `tests/test_splice_binder.py` — Splice helper unit tests (chain rename A->D, residue renumbering, defensive asserts).
- `tests/test_aggregate_mpnn_outputs.py` — MPNN FASTA parsing + MPNN/AF2 seed range + uniqueness tests.
- `tests/test_stage2_designs_halt_gate.py` — Stage 2 designs orchestrator + halt-gate tests (pass-at-boundary + fail paths).
- `tests/test_run_stage2.py` — Stage 2 real-mode helper tests (multimer FASTA writer, target chain sequence extraction, binder-length manifest loader).
- `tests/fixtures/target_3hpj_clean.pdb` — 1-line ATOM stub for stage 0 primary.
- `tests/fixtures/target_2bnr_clean.pdb` — 1-line ATOM stub for stage 0 positive control.
- `tests/fixtures/rfdiffusion/sample.pdb` — Stage 1 cycle-1 stub backbone (legacy).
- `tests/fixtures/stage1/_make_fixtures.py` — Deterministic generator for the Stage 1 mock fixtures.
- `tests/fixtures/stage1/mock_clean.pdb` — Synthetic cleaned PDB with all 8 hotspot Ca residues.
- `tests/fixtures/stage1/mock_target.yaml` — Pre-built TargetManifest for Stage 1 mock smoke tests.
- `tests/fixtures/stage1/mock_design_NNNNN.pdb` — 10 mock backbones (8 pass / 2 fail at threshold 0.50).
- `tests/fixtures/proteinmpnn/sample.fasta` — Stage 2a stub sequence.
- `tests/fixtures/colabfold/sample.pdb` — Stage 2b stub structure.
- `tests/fixtures/colabfold/sample_scores.json` — Stage 2b stub PAE/pLDDT.
- `tests/fixtures/metrics/sample_metrics.json` — Stage 2c stub filter scores.
- `tests/fixtures/stage2/_make_fixtures.py` — Deterministic generator for the 5-control mock PAE/pLDDT JSONs.
- `tests/fixtures/stage2/controls_colabfold/` — Mock ColabFold rank_001 JSON + PDB for P1, P2, P3, N1, N2.
- `tests/fixtures/stage2/stage2_metrics_mock.json` — Mock Stage 2 metrics.json for rule 04 mock branch.
- `tests/fixtures/stage2/designs/_make_fixtures.py` — Generator for the Stage 2 designs mock fixture set (10 backbones, 40 predictions, 4/40 boundary calibration).
- `tests/fixtures/stage2/designs/stage1_backbones/` — 10 mock single-chain Stage 1 binder PDBs.
- `tests/fixtures/stage2/designs/mpnn_outputs/` — 10 mock ProteinMPNN-formatted FASTAs (4 seqs each).
- `tests/fixtures/stage2/designs/af2_predictions/` — 40 ColabFold-shaped prediction subdirs (4 pass / 36 fail at the intermediate cut).
- `tests/fixtures/crosspan/sample_panel.json` — Stage 3 stub.
- `tests/fixtures/embeddings/sample.npz` — Stage 4 stub.
- `tests/fixtures/active_learning/sample_predictions.json` — Stage 5 stub.
- `tests/fixtures/report/sample_manifest.json` — Stage 6 stub.

## notebooks/, reports/

- `notebooks/.gitkeep` — Cycle analysis notebooks land here.
- `reports/.gitkeep` — `cycle_NN.md` reports land here.
