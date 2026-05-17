# CLAUDE.md — pMHC-I Minibinder Design Pipeline

> De novo pMHC-I minibinder design pipeline for CAR-T specificity.
> RFdiffusion → ProteinMPNN → AF2-multimer → cross-pan → ESM-2 diversity → active learning.

This file is the project rulebook. Read it first. For per-stage detail, read
`specs/stageN_*.md`. For a file map, read `INDEX.md`. For personal/local
overrides, read `CLAUDE.local.md` (gitignored).

---

## Locked design decisions

| Decision               | Value                                                  |
|------------------------|--------------------------------------------------------|
| Primary target         | WT1 RMFPNAPYL / HLA-A\*02:01 (PDB **3HPJ**)            |
| Positive control       | Jenkins NY1-B04 vs SLLMWITQC / HLA-A\*02:01 (PDB 2BNR) |
| GPU                    | NVIDIA A100 SXM 80GB (RunPod)                          |
| Orchestrator           | Snakemake 8.x                                          |
| Backbone gen           | RFdiffusion (`rosettacommons/rfdiffusion`)             |
| Sequence design        | ProteinMPNN complex mode, T=0.1, fix MHC chains        |
| Structure prediction   | ColabFold (AF2-multimer, remote MMseqs2 MSA)           |
| Sequence embeddings    | ESM-2 650M, mean-pool over binder residues             |
| Surrogate              | LightGBM + GP baseline                                 |
| Diversity selection    | Farthest Point Sampling (FPS)                          |
| Experiment tracking    | Local JSON + figures committed to repo (no W&B in v0)  |
| Python                 | 3.12+                                                  |
| Dependency manager     | `uv`                                                   |

---

## Two-tab architecture

```
Claude Code Web (authoring)        GitHub (bridge)        RunPod A100 SXM 80GB
  - writes code/specs/configs   →  pmhc_minibinder_   →   /workspace/pipeline/
  - commits + pushes               pipeline               Python 3.12, CUDA 13
  - dry-run gate before push       claude/setup-...       runs real pipeline
```

The user is on Claude Code Web. Actual pipeline execution happens on a RunPod
pod via `git pull` in the pod's Jupyter Lab terminal. **Never assume local
Python/GPU execution.** Code is written for the pod, not for this container.

---

## Workflow rules

1. **Plan mode first.** Each new session enters plan mode, reads the relevant
   `specs/stageN_*.md`, presents a plan, waits for explicit user approval
   ("approved, proceed").
2. **Specs are contracts.** Implement against `specs/stageN_*.md`; never extend
   scope without updating the spec first.
3. **Verify before push.** `snakemake --dry-run --config mock=true -j1` must
   pass. The pre-push hook (`.claude/hooks/pre-push.sh`) enforces this.
4. **Conventional Commits.** `chore:`, `feat:`, `fix:`, `docs:`, `test:`,
   `ci:`, `refactor:`.
5. **`/compact` after each major commit** to preserve session context.
6. **`--mock` everywhere.** Every script supports `--mock` (reads
   `tests/fixtures/`, no GPU, <1s). This is the verification mechanism.
7. **No magic numbers.** All thresholds live in `configs/thresholds.yaml`.
8. **Anchor citations.** Reference papers by anchor name in docstrings:
   `HADRUP_JENKINS_2025`, `BAKER_LAB_2025`, `HOUSEHOLDER_GARCIA_2025`,
   `MARES_IOANNIDIS_2025`, `BENTZEN_HADRUP_2019`.

---

## Filter thresholds (anchored to HADRUP_JENKINS_2025)

| Stage      | Metric         | Cycle 1 (initial) | Final (post partial diffusion) |
|------------|----------------|-------------------|--------------------------------|
| AF2        | iPAE           | < 12 Å            | **< 7 Å**                      |
| AF2        | ipLDDT         | loose             | **> 92**                       |
| AF2        | num_recycles   | 6                 | 6                              |
| ProteinMPNN| NLL            | top 10 %          | top 5 %                        |

Live values are in `configs/thresholds.yaml`. Treat the table above as a
human-readable summary, not as a source of truth.

---

## Controls panel (must run every cycle)

| #  | Control                                    | Expected outcome              |
|----|--------------------------------------------|-------------------------------|
| P1 | Baker WT1 binder vs WT1 / A\*02:01         | iPAE < 7, top 10 % of designs |
| P2 | Jenkins NY1-B04 vs SLLMWITQC / A\*02:01    | iPAE ≈ 6.5 (Fig 1B match)     |
| P3 | Baker WT1 binder vs MART-1 / A\*02:01      | iPAE > P1 + 3 Å               |
| N1 | Scrambled Baker WT1                        | iPAE > 15                     |
| N2 | Random 65-aa (natural AA freqs)            | iPAE > 20                     |

If P1 falls in the bottom 50 % of designs, the pipeline halts and alerts.

---

## Coding conventions

- `pathlib.Path` for all paths; never string concat.
- Type hints everywhere; `mypy --strict` must pass.
- Every script exposes `--mock` and `--config` flags.
- Default to no comments. Comment only when the WHY is non-obvious (cite an
  anchor when the choice comes from a paper).
- Externalize thresholds, paths, hotspot lists — never inline.
- Use pydantic models to load YAML configs.

---

## Pointers

- `INDEX.md` — complete file map (one line per file).
- `specs/stageN_*.md` — stage contracts (Goal, Inputs, Outputs, Tools,
  Anchors, Tasks, Verification, Pitfalls).
- `CLAUDE.local.md` — personal context, RunPod paths (gitignored).
- `bootstrap.sh` — one-command pod setup.
