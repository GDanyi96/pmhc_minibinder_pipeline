# PIPELINE_STATUS.md

Live status of the pMHC-I minibinder pipeline. Updated per cycle / per
architectural decision. For the rulebook see `CLAUDE.md`; for per-stage
contracts see `specs/stageN_*.md`.

## Per-stage status

| Stage | Scope                                | State                  |
|-------|--------------------------------------|------------------------|
| 0     | Target prep (clean PDB, hotspots)    | scaffolded             |
| 1     | RFdiffusion backbones                | **mock-CI green; real-run ready (native driver + weights gated behind `--with-rfdiffusion-weights`)** |
| 2 (controls) | ProteinMPNN + AF2-multimer + P1-N2 panel | mock-CI green; real-run rewired to native install (cycle 1) |
| 2 (designs)  | ProteinMPNN + AF2 fan-in on Stage 1 outputs | **mock-CI green; real-run blocked on pod (post-merge user task)** |
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

## Pointers

- `CLAUDE.md` — rulebook, locked decisions, controls panel summary.
- `INDEX.md` — one-line description per tracked file.
- `bootstrap.sh` — one-command RunPod setup.
- `specs/stageN_*.md` — stage contracts.
