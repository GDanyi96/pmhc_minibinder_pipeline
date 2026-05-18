# PIPELINE_STATUS.md

Live status of the pMHC-I minibinder pipeline. Updated per cycle / per
architectural decision. For the rulebook see `CLAUDE.md`; for per-stage
contracts see `specs/stageN_*.md`.

## Per-stage status

| Stage | Scope                                | State                  |
|-------|--------------------------------------|------------------------|
| 0     | Target prep (clean PDB, hotspots)    | scaffolded             |
| 1     | RFdiffusion backbones                | **mock-CI green; real-run ready (native driver + weights gated behind `--with-rfdiffusion-weights`)** |
| 2     | ProteinMPNN + AF2-multimer + controls | mock-CI green; real-run rewired to native install (cycle 1) |
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

## Pointers

- `CLAUDE.md` — rulebook, locked decisions, controls panel summary.
- `INDEX.md` — one-line description per tracked file.
- `bootstrap.sh` — one-command RunPod setup.
- `specs/stageN_*.md` — stage contracts.
