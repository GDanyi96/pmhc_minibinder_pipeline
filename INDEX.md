# INDEX.md — File map

One-line description per tracked file. New file? Add a line here.

## Top-level

- `CLAUDE.md` — Project rulebook (locked decisions, workflow rules, thresholds).
- `CLAUDE.local.md` — Personal runtime context (gitignored).
- `INDEX.md` — This file.
- `README.md` — Portfolio-facing project overview.
- `LICENSE` — MIT.
- `bootstrap.sh` — One-command RunPod setup (deps, model weights, ProteinMPNN clone, smoke test).
- `PROJECT_STATE.md` — One-screen "where are we": current per-cycle status + non-drifting reference tables (pod, env pins).
- `pyproject.toml` — `uv`-managed Python project + pinned minimal deps.
- `uv.lock` — `uv` resolved dependency lockfile.
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
- `specs/stage3_spec_filter.md` — Placeholder (cycle 04): ProteinMPNN peptide-context specificity filter; distinct from crosspan.
- `specs/stage4_diversity.md` — ESM-2 embeddings + farthest-point sampling diversity selection.
- `specs/stage5_active_learning.md` — LightGBM + GP surrogate, UCB acquisition for next cycle.
- `specs/stage6_reporting.md` — Cycle aggregation, publication figures, decision log.

## workflow/rules/

- `workflow/rules/00_target_prep.smk` — Download/clean target PDBs, emit target.yaml.
- `workflow/rules/01_rfdiffusion.smk` — Run RFdiffusion to generate N backbones (cycle ≤ 02 de-novo path).
- `workflow/rules/01_stage1_subrun_a.smk` — Cycle 03 sub-run A: align BAKER scaffolds + partial diffusion.
- `workflow/rules/01_stage1_subrun_b.smk` — Cycle 03 sub-run B: partial diffusion from the design_2079 hero seed.
- `workflow/rules/01_stage1_merge.smk` — Cycle 03: merge sub-runs A+B into the canonical Stage 1 output + halt gate.
- `workflow/rules/02_proteinmpnn.smk` — Design sequences on each backbone.
- `workflow/rules/02b_proteinmpnn_designs.smk` — Stage 2 designs end-to-end pipeline (splice + MPNN + AF2 fan-in + halt; cycle-02 full 4-chain).
- `workflow/rules/02c_stage2_subrun_a.smk` — Cycle 03 Stage 2 sub-run A: truncated 3-chain designs (`--target-layout truncated`); explicit-target rule, not in `rule all`.
- `workflow/rules/03_colabfold.smk` — AF2-multimer prediction via ColabFold.
- `workflow/rules/03b_af2_designs.smk` — Marker-only DAG anchor downstream of stage2_designs.
- `workflow/rules/04_metrics.smk` — Compute iPAE, ipLDDT, NLL filters.
- `workflow/rules/05_crosspan.smk` — Cross-reactivity panel against off-targets.
- `workflow/rules/06_embedding.smk` — ESM-2 embeddings + FPS diversity selection.
- `workflow/rules/07_active_learning.smk` — Surrogate training + UCB acquisition.
- `workflow/rules/08_reporting.smk` — Aggregate to cycle report.

## workflow/scripts/

- `workflow/__init__.py` — Marks workflow as a package.
- `workflow/scripts/__init__.py` — Marks workflow.scripts as a package (mypy + imports).
- `workflow/scripts/prep_target.py` — Stage 0 entry point.
- `workflow/scripts/prep_baker_target.py` — Cycle 03 sub-run A: truncate 3hpj_clean.pdb to the BAKER-format groove-only target (chain B=HLA[1:180], C=peptide; no β2m/α3).
- `workflow/scripts/run_rfdiffusion.py` — Stage 1 entry point.
- `workflow/scripts/align_scaffolds.py` — Cycle 03: BioPython superposition of BAKER scaffolds onto the reference (align_chainB.py equivalent); `baker_layout=True` dispatches to align_baker_scaffolds.
- `workflow/scripts/align_baker_scaffolds.py` — Cycle 03 sub-run A: align BAKER fused-chain scaffolds onto the truncated target, emit A=binder/B=HLA/C=peptide.
- `workflow/scripts/partial_diffuse.py` — Cycle 03: RFdiffusion partial-diffusion wrapper (mock synthesizes geometry-passing 4-chain designs).
- `workflow/scripts/contact_filter.py` — Cycle 03: BioPython binder↔peptide C-beta contact gate.
- `workflow/scripts/setup_cycle03_inputs.py` — Pod-only one-shot: symlink BAKER scaffolds + stitch the design_2079 seed complex.
- `workflow/scripts/run_proteinmpnn.py` — Stage 2a: ProteinMPNN complex-mode wrapper. **DEPRECATED for real-mode** (cycle 02+); mock paths retained for cycle-1 controls. Real-mode bypasses this and invokes upstream ProteinMPNN directly from `scripts/run_stage2.py`.
- `workflow/scripts/run_colabfold.py` — Stage 2b: ColabFold (AF2-multimer) native wrapper.
- `workflow/scripts/compute_metrics.py` — Stage 2c: iPAE / ipLDDT / BSA computation (interface-8 Å, decomposed pep/mhc; truncated-layout aware).
- `workflow/scripts/splice_binder.py` — Stage 2 designs: splice Stage 1 binder onto cleaned pMHC (cycle-02 A→D; truncated sub-run-A fork A=HLA/B=peptide/C=binder).
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
- `scripts/run_stage2.py` — Stage 2 designs orchestrator (splice + MPNN + AF2 fan-in + halt gate; `--target-layout` full/truncated).
- `scripts/run_stage1_subrun.py` — Cycle 03: drives one partial-diffusion sub-run (A or B); emits per-subrun designs.jsonl + summary.
- `scripts/merge_stage1_subruns.py` — Cycle 03: merge sub-runs A+B into canonical stage1 output + geometry halt gate.

## analysis/scripts/

- `analysis/scripts/12_design3010_peptide_contacts.py` — Cycle 03 hero (design_3010_seq00): per-peptide-residue min-distance + closest binder residue + PAE (truncated 3-chain layout).
- `analysis/scripts/13_design2079_peptide_contacts.py` — Cycle 02 hero (design_2079_seq00): same per-residue peptide-contact analysis on the full 4-chain layout (confirms peptide-blind framework binder).

## configs/

- `configs/target_wt1_a0201.yaml` — Locked primary target spec (chains, hotspots, length range 70-110).
- `configs/target_2bnr_a0201.yaml` — Positive-control target spec (Jenkins NY1-B04 / 2BNR).
- `configs/thresholds.yaml` — All numerical filter thresholds (iPAE, ipLDDT, NLL, etc).
- `configs/rfdiffusion_default.yaml` — RFdiffusion noise scaling, num_designs, partial_T (contigmap built at runtime).
- `configs/seeds.yaml` — Stage 1 seed formula (`cycle * 1000 + design_index`) + per-cycle reserved ranges.
- `configs/proteinmpnn_chains.json` — Per-PDB chain assignment for ProteinMPNN complex mode.
- `configs/proteinmpnn_default.yaml` — Stage 2 designs MPNN sampling defaults (T=0.1, 4 seqs/backbone, chain D designed).
- `configs/af2_stage2.yaml` — Stage 2 designs AF2 fan-in config (cycle 03: top-100, num_recycles=6, intermediate cut iPAE<10/ipLDDT>88, halt 0.10).
- `configs/rfdiffusion_subrun_a.yaml` — Cycle 03 sub-run A partial-diffusion config (BAKER scaffolds, partial_T=15).
- `configs/rfdiffusion_subrun_b.yaml` — Cycle 03 sub-run B partial-diffusion config (design_2079 seed, partial_T=10).
- `configs/proteinmpnn_cycle03.yaml` — Cycle 03 MPNN config; adds bias_AA_jsonl to counter the Ala-heavy prior (Trap #30).
- `configs/proteinmpnn_bias_aa.json` — ProteinMPNN --bias_AA_jsonl: A:-2.0, E/L/R:+1.0.

## data/

- `data/targets/.gitkeep` — PDB files land here (downloaded by bootstrap.sh).
- `data/controls/controls_manifest.yaml` — P1, P2, P3, N1, N2 control definitions.
- `data/scaffolds/baker_library/.gitkeep` — BAKER scaffolds symlinked here pod-side (setup_cycle03_inputs.py); contents gitignored.
- `data/seeds/.gitkeep` — design_2079 seed complex stitched here pod-side; contents gitignored.

## docs/

- `docs/narrative.md` — Portfolio scientific narrative (motivation, target, funnel, cross-cycle framing, conclusions).
- `docs/cycle_02.md` — Cycle 02 detailed results: de novo RFdiffusion campaign (full 4-chain target).
- `docs/cycle_03.md` — Cycle 03 detailed results: BAKER-scaffold partial-diffusion campaign (truncated target, sub-run A).
- `docs/methodological_lessons.md` — Science-bearing conceptual lessons (peptide-vs-framework blindness, iPAE definition audit, charge hypothesis).
- `docs/pod_quickstart.md` — Fresh-pod three-line sequence, debug flags, recovery from failed metrics steps.
- `docs/known_traps.md` — Empirical engineering gotchas bringing the pipeline up on RunPod (Docker, uv, ColabFold layout, PAE key, JAX drift, LD_LIBRARY_PATH, AF2 specificity blind spot).
- `docs/figures/.gitkeep` — Committed publication figures land here.
- `docs/figures/cycle02_funnel.png` — Cycle 02 design funnel.
- `docs/figures/cycle02_ipae_distribution.png` — Cycle 02 iPAE distribution.
- `docs/figures/cycle02_stage1_geometry.png` — Cycle 02 Stage 1 geometry-pass distribution.
- `docs/figures/cycle03_funnel.png` — Cycle 03 design funnel.
- `docs/figures/design2079_aa_composition.png` — Cycle 02 hero amino-acid composition (Ala-heavy prior).
- `docs/figures/design2079_lottery.png` — design_2079 sequence-lottery analysis.
- `docs/figures/design_3010_peptide_interface.png` — Cycle 03 hero design_3010 peptide interface.
- `docs/figures/pep_vs_framework_scatter.png` — `iface_pep` vs `iface_mhc` scatter (specificity axis).

## results/

- `results/master_design_journey.csv` — Cross-cycle per-design master table (decomposed iPAE, MPNN, composition, scaffold lineage, peptide-contact tier).
- `results/cycle_01/stage2/metrics.json` — Cycle 01 Stage 2 control metrics.
- `results/cycle_01/stage2/metrics_A100_baseline.json` — Cycle 01 Stage 2 control metrics (A100 baseline run).
- `results/cycle_01/stage2/real_run.log` — Cycle 01 ColabFold real-run log.
- `results/cycle_02/stage1/.gitkeep` — Cycle 02 Stage 1 results dir placeholder.
- `results/cycle_02/stage2/.gitkeep` — Cycle 02 Stage 2 results dir placeholder.
- `results/cycle_03/analysis/candidate_dossier.csv` — Cycle 03 ranked candidate dossier (peptide-contact tier + interface metrics + scaffold/composition).
- `results/cycle_03/analysis/scaffold_lineage.csv` — Cycle 03 BAKER scaffold lineage (HLA Cα RMSD, cluster).
- `results/cycle_03/analysis/design_3010_seq00_af2_model.pdb` — Cycle 03 hero AF2 model.

## tests/

- `tests/__init__.py` — Marks tests as a package.
- `tests/test_dry_run.py` — Asserts `snakemake --dry-run --config mock=true -j1` exits 0.
- `tests/test_halt_gate.py` — Stage 2 controls halt-gate unit tests.
- `tests/test_compute_metrics.py` — Stage 2 metrics computation tests.
- `tests/test_run_rfdiffusion.py` — Stage 1 helper unit tests (contigmap, seeds, geometry pass).
- `tests/test_run_rfdiffusion_real_writer.py` — Regression: real RFdiffusion `design_{seed}.pdb` writer/reader filename contract (cycle 02 n_completed=0 incident).
- `tests/test_stage1_halt_gate.py` — Stage 1 orchestrator + halt-gate tests (pass + fail paths).
- `tests/test_splice_binder.py` — Splice helper unit tests (chain rename A->D, residue renumbering, defensive asserts).
- `tests/test_aggregate_mpnn_outputs.py` — MPNN FASTA parsing + MPNN/AF2 seed range + uniqueness tests.
- `tests/test_stage2_designs_halt_gate.py` — Stage 2 designs orchestrator + halt-gate tests (pass-at-boundary + fail paths).
- `tests/test_stage2_subrun_a.py` — Cycle 03 Stage 2 sub-run A truncated path (3-chain splice + MPNN chain C + truncated metrics).
- `tests/test_align_scaffolds.py` — Cycle 03 align_scaffolds superposition tests.
- `tests/test_prep_baker_target.py` — Cycle 03 sub-run A: truncation CA-count / layout / numbering tests.
- `tests/test_align_baker_scaffolds.py` — Cycle 03 sub-run A: BAKER fused-chain alignment + dispatch tests.
- `tests/test_controls_truncated.py` — Cycle 03: truncated controls FASTA + target-layout iPAE decomposition tests.
- `tests/test_partial_diffuse.py` — Cycle 03 partial_diffuse mock synthesis + geometry-pass tests.
- `tests/test_contact_filter.py` — Cycle 03 contact_filter peptide-contact gate tests.
- `tests/test_af2_metrics_decomposition.py` — Cycle 03 decomposed iPAE (peptide/MHC) sub-metric tests.
- `tests/test_stage1_subruns.py` — Cycle 03 sub-run + merge orchestrator tests (pass + fail halt paths).
- `tests/test_run_stage2.py` — Stage 2 real-mode helper tests (multimer FASTA writer, target chain sequence extraction, binder-length manifest loader).
- `tests/test_no_absolute_paths_in_committed_fixtures.py` — Guard: no committed fixture YAML holds an absolute (`/`-prefixed) path.

### tests/fixtures/

- `tests/fixtures/target_3hpj_clean.pdb` — 1-line ATOM stub for stage 0 primary.
- `tests/fixtures/target_2bnr_clean.pdb` — 1-line ATOM stub for stage 0 positive control.
- `tests/fixtures/targets/3hpj_baker_truncated_mock.pdb` — Mock BAKER-format truncated target (chain B=HLA[180], C=peptide[9]).
- `tests/fixtures/rfdiffusion/sample.pdb` — Stage 1 cycle-1 stub backbone (legacy).
- `tests/fixtures/baker_library_mock/_make_fixtures.py` — Generator for the 3 mock BAKER scaffolds + the stitched design_2079 seed.
- `tests/fixtures/baker_library_mock/scaf{0,1,2}.pdb` — 3 synthetic mini BAKER scaffolds (binder chain A + truncated target chain B).
- `tests/fixtures/baker_truncated_mock/_make_fixtures.py` — Generator for the BAKER-truncated mock target + fused-chain mock scaffolds.
- `tests/fixtures/baker_truncated_mock/scaf{0,1}.pdb` — Mock BAKER scaffolds (binder chain A + fused chain B = HLA[180]+peptide[9]).
- `tests/fixtures/design_2079_mock_seed.pdb` — Synthetic 4-chain design_2079 seed complex (A/B/C motif + chain D binder).
- `tests/fixtures/stage1/.gitkeep` — Stage 1 fixtures dir placeholder.
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
- `tests/fixtures/stage2/designs/.gitkeep` — Stage 2 designs fixtures dir placeholder.
- `tests/fixtures/stage2/designs/_make_fixtures.py` — Generator for the Stage 2 designs mock fixture set (10 backbones, 40 predictions, 4/40 boundary calibration).
- `tests/fixtures/stage2/designs/mock_pmhc.pdb` — Mock cleaned full 4-chain pMHC target.
- `tests/fixtures/stage2/designs/mock_target.yaml` — TargetManifest pointing at mock_pmhc.pdb.
- `tests/fixtures/stage2/designs/stage1_backbones/` — 10 mock single-chain Stage 1 binder PDBs.
- `tests/fixtures/stage2/designs/stage1_designs.jsonl` — Mock Stage 1 designs manifest.
- `tests/fixtures/stage2/designs/stage1_summary.json` — Mock Stage 1 summary.
- `tests/fixtures/stage2/designs/mpnn_outputs/` — 10 mock ProteinMPNN-formatted FASTAs (4 seqs each).
- `tests/fixtures/stage2/designs/af2_predictions/` — 40 ColabFold-shaped prediction subdirs (4 pass / 36 fail at the intermediate cut).
- `tests/fixtures/stage2/designs_truncated/_make_fixtures.py` — Generator for the Stage 2 sub-run A (truncated 3-chain) mock fixture set.
- `tests/fixtures/stage2/designs_truncated/mock_truncated_pmhc.pdb` — Mock cleaned BAKER truncated target (B=HLA, C=peptide).
- `tests/fixtures/stage2/designs_truncated/mock_truncated_target.yaml` — TargetManifest pointing at mock_truncated_pmhc.pdb.
- `tests/fixtures/stage2/designs_truncated/stage1_backbones/` — 3 mock 3-chain Stage 1 backbones (A=binder, B=HLA, C=peptide).
- `tests/fixtures/stage2/designs_truncated/stage1_designs.jsonl` — Mock truncated Stage 1 designs manifest.
- `tests/fixtures/stage2/designs_truncated/subrun_summary.json` — Mock sub-run A summary.
- `tests/fixtures/stage2/designs_truncated/mpnn_outputs/` — 3 mock ProteinMPNN FASTAs (chain C designed).
- `tests/fixtures/stage2/designs_truncated/af2_predictions/` — Mock ColabFold prediction subdirs in HLA:peptide:binder order.
- `tests/fixtures/crosspan/sample_panel.json` — Stage 3 stub.
- `tests/fixtures/embeddings/sample.npz` — Stage 4 stub.
- `tests/fixtures/active_learning/sample_predictions.json` — Stage 5 stub.
- `tests/fixtures/report/sample_manifest.json` — Stage 6 stub.
