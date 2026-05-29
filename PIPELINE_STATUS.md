# PIPELINE_STATUS.md

Live status of the pMHC-I minibinder pipeline. Updated per cycle / per
architectural decision. For the rulebook see `CLAUDE.md`; for per-stage
contracts see `specs/stageN_*.md`.

## Per-stage status

| Stage | Scope                                | State                  |
|-------|--------------------------------------|------------------------|
| 0     | Target prep (clean PDB, hotspots)    | scaffolded             |
| 1     | RFdiffusion backbones                | ✓ **real runs complete — cycle 02 de novo (200 designs, 13% pass) + cycle 03 partial diffusion sub-run A (150 designs, 72/150 = 48% pass)** |
| 2 (controls) | ProteinMPNN + AF2-multimer + P1-N2 panel | ✓ validated — full target (cycle 1) + truncated target (cycle 3) |
| 2 (designs)  | ProteinMPNN + AF2 fan-in on Stage 1 outputs | ✓ **real runs complete — cycle 02 full target (50 folds, hero design_2079) + cycle 03 truncated sub-run A (288 folds; 1 peptide reader, design_3010)** |
| 3     | Cross-pan off-target grid            | not started            |
| 4     | ESM-2 embeddings + FPS               | not started            |
| 5     | LightGBM + GP active learning        | not started            |
| 6     | Cycle reporting                      | not started            |

## CI status

- Mock-mode CI is green on `claude/stage-2-validation-pipeline-uIoLz` as of
  2026-05-17: `ruff`, `black --check`, `mypy --strict`, `snakemake --dry-run
  --config mock=true -j1`, and `pytest -q` all pass.
- CI installs via `uv sync --all-extras` (one source of truth =
  `pyproject.toml`).

## Real-run status

- **Stage 2** real run was rewired to native install in cycle 1
  (`run_colabfold.py` and `run_proteinmpnn.py` shell out natively; no
  Docker).
- **Stage 1** (RFdiffusion) is now native too: `workflow/scripts/run_rfdiffusion.py`
  shells out to `/workspace/RFdiffusion/scripts/run_inference.py` via a
  module-level `_SEED_THREADING_MODE` constant (defaulting to the safe
  per-design subprocess pattern; the operator flips to
  `single_subprocess` after the zero-LOC pod recon documented in the
  module docstring). Weights are gated behind
  `bootstrap.sh --with-rfdiffusion-weights` with sha256 verification.

## Locked architectural decisions

| Decision                                                                    | Rationale |
|-----------------------------------------------------------------------------|-----------|
| **Native install on RunPod pod; the pod IS the container, no DinD**         | RunPod pods are containers themselves; Docker-in-Docker is not configured and won't be. All heavy tools (ColabFold, ProteinMPNN, RFdiffusion) install natively into the pod's Python env or as cloned upstream repos. |
| ColabFold via the `[colabfold]` optional extra (sokrypton/ColabFold git URL) | Heavy (~5 GB) — production `uv sync` should not pull it. Pinning a commit hash is queued for after a known-good pod revision is validated. |
| ProteinMPNN via upstream clone at `/workspace/ProteinMPNN`                  | No pip wheel exists for `dauparas/ProteinMPNN`. `bootstrap.sh` clones it; `run_proteinmpnn.py` invokes `protein_mpnn_run.py` via `sys.executable`. Path overridable via `PROTEINMPNN_DIR`. |
| Dev tools (`ruff`, `black`, `mypy`, `pytest`, `types-PyYAML`) live in `[dev]` | Production `uv sync` stays lean. `bootstrap.sh` and CI use `uv sync --all-extras`. |

## Cycle 1 — Stage 2 — H100 recalibration (2026-05-18)

2026-05-18: Cycle 1 controls re-run on H100 SXM US-NE-1 (pod brief_beige_mole)
after the unplanned cross-region migration from US-CA-2. All five controls
pass; halt gate verdict PASS. Positives drift ≤ +0.34 Å iPAE vs A100
baseline; ipLDDT shifts ≤ 4.4 points; BSA shifts ≤ 664 Å² (N2). metrics.json
now reflects H100 numbers (canonical). Original A100 numbers preserved as
metrics_A100_baseline.json sidecar in the same directory. Note: the baseline
file preserves the historical false-halt verdict from the old rank-based halt
rule; the underlying numbers passed under the current biology-correct rule.
Halt gate margin intact: iPAE gap 20.0 Å (threshold ≥10), ipLDDT gap 59.8
(threshold ≥30). See trap #20 for the GPU-drift envelope.

## Cycle 03 prep (2026-05-21)

Prep-only PR (`claude/cycle-03-prep-g2Eu9`) wiring the cycle-03 architecture;
no GPU / no real run yet. Addresses cycle-02's placement deficit (65% of
designs made zero hotspot contact, one hero out of 200) by switching Stage 1
from de-novo RFdiffusion to **partial diffusion** from two seeded sub-runs.

- **Stage 1 mode switch**: `Snakefile` now derives `stage1_mode` from the
  cycle (`partial` when cycle >= 3, else `denovo`; overridable via
  `--config stage1_mode=`). A conditional `ruleorder` selects either the
  legacy `rule rfdiffusion` or the new `stage1_merge` (sub-runs A+B), both of
  which emit the canonical `stage1/rfdiffusion/{stage1_summary.json,
  designs.jsonl}`. Cycle 01/02 paths are unchanged.
- **Sub-run A** (`rule stage1_subrun_a`): BAKER scaffold library →
  `align_scaffolds.py` (BioPython superposition onto `3hpj_clean.pdb`, a
  Rosetta-free `align_chainB.py` equivalent) → `partial_diffuse.py`
  (`partial_T=15`). **Sub-run B**: cycle-02 hero seed complex →
  `partial_diffuse.py` (`partial_T=10`). `partial_T` values are conservative
  blind defaults (recalibrate in cycle 04).
- **Stage 2 cycle-03 changes**: Ala-counter-bias MPNN config
  (`proteinmpnn_cycle03.yaml` + `proteinmpnn_bias_aa.json`, Trap #30); a
  BioPython peptide-contact gate (`contact_filter.py`, C-beta ≤ 5 Å); AF2
  iPAE decomposed into `ppi_pae_int_peptide` / `ppi_pae_int_mhc`
  (`compute_metrics.py`); halt tightened to iPAE ≤ 10, `num_recycles=6`,
  `fan_in_top_n=100` (`af2_stage2.yaml`).
- **Pod-only artefacts** (BAKER scaffolds, `design_2079` hero) stay gitignored
  and are materialized on the pod by `setup_cycle03_inputs.py` (symlink
  scaffolds + stitch the hero chain D onto `3hpj_clean.pdb`). CC develops
  against mock fixtures (`tests/fixtures/baker_library_mock/`,
  `tests/fixtures/design_2079_mock_seed.pdb`).
- **Verification**: `snakemake --config mock=true cycle=99 -j1` runs the full
  partial-diffusion DAG end to end (sub-runs → merge → Stage 2 → report),
  `cycle=01` legacy path stays green, ruff/black/mypy/pytest all pass.
- **Deferred to cycle 04**: the ProteinMPNN peptide-context specificity
  (`mpnn_spec_filter`) filter — see `specs/stage3_spec_filter.md`.

## Cycle 03 BAKER truncation (sub-run A prep)

- **Why**: BAKER scaffolds carry the target as a fused chain B = HLA[1:180] +
  peptide (binder on A). Aligning that onto `3hpj_clean.pdb` (chain B = β2m)
  silently superposes HLA onto β2m → garbage geometry, ~6 h A100 wasted
  (Trap #31). Sub-run A needs a target whose chain layout matches BAKER's
  groove-only fragment.
- **Truncated target**: `prep_baker_target.py` →
  `data/targets/3hpj_baker_truncated.pdb` (chain B = HLA α1+α2 res 1-180,
  chain C = peptide 1-9; β2m and α3 dropped; CA counts 180 B / 9 C). Built by
  `rule prep_baker_target` (mock copies the fixture) and pod-side by
  `setup_cycle03_inputs.py`. **Sub-run B keeps the full target** (`design_2079`
  was AF2-validated against 4-chain `3hpj_clean.pdb`).
- **Alignment**: `align_baker_scaffolds.py` superposes the scaffold's chain-B
  HLA substring onto the truncated reference and rewrites to our
  `A=binder / B=HLA / C=peptide` layout; `align_scaffolds(..., baker_layout=True)`
  dispatches to it. `rule stage1_subrun_a` passes the truncated target as the
  align `--reference-pdb`; the geometry/motif reference stays the full manifest.
  The mock DAG keeps the lightweight aligner + `mock_clean.pdb`, so cycle 99/01
  stay green.
- **Controls recalibration**: `run_controls.py --target {full,truncated}`. The
  truncated baseline folds binders against the β2m-free groove
  (`HLA : peptide : binder` FASTA → 3-chain A=HLA, B=peptide, C=binder) and
  writes `results/cycle_NN/controls_truncated_baseline/metrics_truncated_baseline.json`.
  Launch gate (enforced when run on pod): P1 truncated iPAE within 1.0 Å of
  cycle-01 P1 (4.5 ± 1.0) and N1 truncated iPAE > 15, else ABORT sub-run A.
- **Metrics adapter**: `compute_metrics.LAYOUT_CHAINS[target_layout]` selects
  binder/peptide/MHC chains (full: D / (C) / (A,B); truncated: C / (B) / (A)),
  threaded through `metrics_for_design` / `metrics_for_control` and
  `run_stage2._compute_metrics_for_predictions`. Default `full` → no behavior
  change for existing callers.
- **Verification**: `prep_baker_target` mock + 180/9 unit tests; align mock
  layout test; `snakemake --dry-run cycle=03` and `cycle=99` e2e green; cycle
  01/02 legacy green; ruff/black/mypy/pytest all pass.

## Cycle 03 — executed + analysis (2026-05-29)

The two sections above are **prep**. Sub-run A then ran end-to-end on the A100 pod and is **merged to `main`**.

- **Stage 1**: 152 BAKER scaffolds aligned → 150 partial-diffused (`partial_T=15`) → **72/150 (48%) geometry-pass** (cycle 02 was 13%). Scaffold HLA-CA RMSD bimodal (56 native < 1 Å, ~95 cross-allele 4–6 Å); cross-allele scaffolds passed *more* often (54% vs 38%), so HLA-RMSD is **not** predictive of transfer (Trap #40).
- **Stage 2**: 288 AF2 folds (72 × 4 MPNN seqs, `num_recycles=6`, AA bias A:−2.0 / E,L,R:+1.0). **Only 91/288 (32%) engage the peptide.** One control-grade reader: **`design_3010_seq00`** (`iface_pep` 2.10) reads WT1 N5 + Y8 through binder R55. The best overall interface (`design_3084_seq02`, ipTM 0.89) and the framework champions are *not* peptide readers.
- **Metric audit**: cycle 02 stored iPAE as interface-8 Å (Jenkins Fig 1B), cycle 03 as the Bennett position-slice; both recomputed. Canonical ranking is now **interface-8 Å**, decomposed into `iface_pep` / `iface_mhc` (the specificity axis). See `CLAUDE.md` → *Metric conventions*.
- **Cycle-02 hero reclassified**: `design_2079_seq00` is **peptide-blind** (`iface_pep` ∞; binder 28–40 Å from the peptide) — a confident MHC-framework binder, not a peptide reader.
- **Diagnosis**: framework bias is **structural** (RFdiffusion placement), not compositional — a charge-driven hypothesis was tested and falsified.
- **Committed**: `results/master_design_journey.csv`, `results/cycle_03/analysis/`, `analysis/scripts/{12,13}_*.py`, `docs/figures/`, and the portfolio docs (`docs/narrative.md`, `docs/cycle_02.md`, `docs/cycle_03.md`, `docs/methodological_lessons.md`).

## Cycle 04 — not started

Scope TBD (not detailed here yet). Carries the structural-placement diagnosis and the `iface_pep` specificity gate as the ranking convention into the next round.

## Pointers

- `CLAUDE.md` — rulebook, locked decisions, controls panel summary.
- `INDEX.md` — one-line description per tracked file.
- `bootstrap.sh` — one-command RunPod setup.
- `specs/stageN_*.md` — stage contracts.
