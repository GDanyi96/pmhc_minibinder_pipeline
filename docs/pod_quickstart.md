# RunPod quickstart

The pod IS the container (no Docker-in-Docker). All heavy tools install
natively into the pod's Python env via `uv` + a `dauparas/ProteinMPNN`
clone. See `docs/known_traps.md` for the empirical gotchas this layout
avoided/encountered in cycle 1.

## Fresh-pod three-line sequence

```bash
cd /workspace/pipeline
git pull
bash bootstrap.sh
uv run python scripts/run_controls.py        # ~10 min on A100, real ColabFold
```

`bootstrap.sh` is idempotent — every step prints `… already present,
skipping` when its outputs exist, so reruns are cheap. With
`--with-rfdiffusion-weights` it also fetches Stage 1 weights (defer until
Stage 1).

## Debug loop (no GPU cost)

Iterating on halt rules, metric definitions, or controls schema doesn't
need a fresh ColabFold pass. ColabFold's own "outputs already exist?"
idempotency check still costs ~35 s/iteration on a 5-control set. Skip
it explicitly:

```bash
uv run python scripts/run_controls.py --metrics-only        # rerun metrics+halt only
uv run python scripts/run_controls.py --skip-colabfold      # rerun FASTA+metrics+halt
```

`--metrics-only` skips: negatives generation, FASTA writing, and the
ColabFold call. It only re-reads the existing `colabfold/controls/{P1..N2}/`
trees and recomputes metrics + halt gate.

`--skip-colabfold` is a softer variant: still regenerates FASTAs and N1/N2,
just skips the ColabFold subprocess.

## Mock smoke test (no GPU, no network)

```bash
uv run python scripts/run_controls.py --mock --cycle 99
```

Copies pre-baked ColabFold fixtures into `results/cycle_99/`. Halt gate
should pass; ipLDDT / iPAE values come from
`tests/fixtures/stage2/controls_colabfold/`. Use a non-`01` cycle if you
don't want to clobber committed cycle-1 evidence.

## `LD_LIBRARY_PATH` warning

RunPod's base image sets `LD_LIBRARY_PATH=/usr/local/cuda/lib64` (CUDA 13).
JAX cu12 wheels bundle CUDA 12 libs; the system path override SIGSEGVs
cuDNN. `bootstrap.sh` unsets it for the bootstrap shell and all children;
`scripts/run_controls.py` also pops it at import time. **But**: a new shell
opened after bootstrap completes (e.g. another Jupyter terminal, an
`ssh` re-login) will re-pick the variable from `/etc/environment`. Before
running `colabfold_batch` or any JAX-importing script in such a shell,
`unset LD_LIBRARY_PATH` first.

## Recovery from a failed metrics step

If the controls run reaches ColabFold successfully but fails during
metrics extraction or the halt gate, the ColabFold outputs are on disk at:

```
results/cycle_01/stage2/colabfold/controls/{P1,P2,P3,N1,N2}/
```

`workflow/scripts/run_colabfold.py` symlinks ColabFold 1.6.1's flat
outputs into these per-id directories after the subprocess returns, so
`compute_metrics.py` finds them via its existing glob. To recover:

```bash
# Fix whatever broke, then:
uv run python scripts/run_controls.py --skip-colabfold
```

No GPU time wasted, no rebuilt FASTAs, ~5 s total.
