# PIPELINE_STATUS.md

Live status of the pMHC-I minibinder pipeline. Updated per cycle / per
architectural decision. For the rulebook see `CLAUDE.md`; for per-stage
contracts see `specs/stageN_*.md`.

## Per-stage status

| Stage | Scope                                | State                  |
|-------|--------------------------------------|------------------------|
| 0     | Target prep (clean PDB, hotspots)    | scaffolded             |
| 1     | RFdiffusion backbones                | stub (real run blocked on Stage 1 native install) |
| 2     | ProteinMPNN + AF2-multimer + controls | **mock-CI green; real-run rewired to native install (this PR)** |
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

- **Stage 2** real run on the RunPod A100 pod previously failed because
  `run_colabfold.py` and `run_proteinmpnn.py` shelled out to
  `docker run …`, and the pod has no Docker (`bash: docker: command not
  found`). This PR replaces both with native subprocess calls.
- Stage 1 (RFdiffusion) remains Docker-based in code; rewiring is queued as
  a follow-up. Spec text is already aligned in this PR.

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
