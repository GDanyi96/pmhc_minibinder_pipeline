# Known traps

Gotchas encountered while bringing the pipeline up on RunPod, with the
fix or workaround for each. Add an entry when something bites for the
second time.

## Docker on RunPod (the pod IS a container)

`bash: docker: command not found`. RunPod pods are themselves containers;
Docker-in-Docker is not configured. All heavy tools must install
natively into the pod's Python env (`[colabfold]` extra,
`/workspace/ProteinMPNN` clone). The original `docker run colabfold/…`
and `docker run proteinmpnn:latest …` invocations were dead on the pod.
**Fix**: native install architecture, locked in `CLAUDE.md`.

## System Python vs `uv` venv

The pod's system Python (`/usr/bin/python3.12`) has none of our deps. Bare
`python scripts/run_controls.py` will `ImportError` on `Bio.PDB`, yaml,
pydantic, etc. **Always use `uv run`** — it resolves the project's `.venv`
created by `uv sync`.

## `uv sync` evicts dev tools

Default `uv sync` installs only `[project.dependencies]`. ruff, black,
mypy, pytest, types-PyYAML live in `[project.optional-dependencies.dev]`
and disappear from `.venv` on a plain `uv sync`. **Fix**:
`uv sync --all-extras` in `bootstrap.sh` (installs dev + colabfold +
proteinmpnn); `uv sync --extra dev` in CI (skip the heavy git URLs CI
doesn't need).

## `git checkout` aborts on untracked `uv.lock`

If `uv.lock` is untracked locally but committed on the branch you're
switching to, `git checkout` refuses to overwrite it. **Fix**: commit
`uv.lock` (the project policy per `.gitignore`'s `# uv.lock` comment).

## ColabFold 1.6.1 writes outputs flat

CF 1.6.1 writes `{out_dir}/{id}_*.{json,pdb,a3m}`, **not** per-id
subdirectories. `compute_metrics.find_artifact` globs
`{out_dir}/{id}/...` and finds nothing → `FileNotFoundError`. **Fix**:
`workflow/scripts/run_colabfold.py:_reshape_flat_outputs` symlinks each
id's flat files into `{out_dir}/{id}/` after the subprocess returns.
Original files stay in place for debugging.

## PAE key rename (`predicted_aligned_error` → `pae`)

ColabFold ≥ 1.6.0 emits the per-residue PAE matrix under the key `pae`;
older versions used `predicted_aligned_error`. **Fix**: `load_pae` tries
`pae` first, falls back to `predicted_aligned_error` for legacy fixtures
and pre-1.6 installs. Mock fixture generator emits `pae` going forward.

## JAX version drift breaks ColabFold

ColabFold 1.6.1 bundles alphafold that calls `jnp.clip(a_max=…)` from
`modules_multimer.py:_relative_encoding`. The `a_max` kwarg was removed
in JAX 0.5+; unpinned uv resolution lands on 0.10.x and the
`colabfold_batch` subprocess raises:

```
TypeError: clip() got an unexpected keyword argument 'a_max'
```

**Fix**: pin `jax[cuda12]==0.4.34` + `jaxlib==0.4.34` in the
`[colabfold]` extra in `pyproject.toml`.

## `LD_LIBRARY_PATH` override forces wrong CUDA libs

RunPod's base image sets `LD_LIBRARY_PATH=/usr/local/cuda/lib64`
(CUDA 13). JAX cu12 wheels bundle their own CUDA 12 libs; the system
path override loads mismatched libs into the JAX process →
`SIGSEGV` in cuDNN at the first forward pass. **Fix**: `unset
LD_LIBRARY_PATH` at the top of `bootstrap.sh`, and
`os.environ.pop("LD_LIBRARY_PATH", None)` at the top of
`scripts/run_controls.py`. The colabfold subprocess inherits the cleaned
env. A new shell opened post-bootstrap will re-pick the variable from
`/etc/environment`; unset manually before running JAX-importing scripts
from such a shell.

## AF2-multimer specificity blind spot

AF2-multimer over-predicts pMHC binding to non-cognate targets for
designed peptide-targeting binders. Confirmed empirically in cycle 1
controls — P3 (wt1-5 vs MART-1) scored iPAE 4.54 Å, essentially
indistinguishable from P1 (wt1-5 vs cognate WT1, iPAE 4.88 Å). Stage 2
iPAE cannot assess specificity; this is a known limitation of
AF2-multimer for designed binders, not a pipeline bug
(HOUSEHOLDER_GARCIA_2025, MARES_IOANNIDIS_2025).

**Fix**: removed the P3-based halt rule. Specificity is assessed in
Stage 3 cross-panning against the off-target peptide grid + Hamming
proteome scan, not in Stage 2. The `enforce_halt_gate` function still
logs a single INFO-level line comparing P3 vs P1 so the cognate-vs-
non-cognate gap is recorded in the run output, but it does not affect
the halt decision.

## Stage 1 geometry-pass threshold is calibration-only for cycle 2

Trap #13. `scripts.run_stage1.HALT_THRESHOLD = 0.50` (in
`configs/seeds.yaml`'s neighborhood; per-stage halt rule lives in
`scripts/run_stage1.py`). 0.50 is a "not catastrophically broken" floor,
not a published threshold. We have zero empirical data on what fraction
of RFdiffusion outputs pass our specific geometry checks (>=3 hotspot
Ca contacts within 10 A + length in [length_min-5, length_max+5] + no
internal Ca-Ca clash <3.5 A) for our specific contigmap. **Fix**: log
the empirical `fraction_geometry_pass` in `stage1_summary.json` on every
cycle and recalibrate the threshold in cycle 3 once the cycle-2
distribution is observed. Tighten toward HADRUP_JENKINS_2025's de novo
yield bracket if cycle-2 data supports it. Documented in
`specs/stage1_rfdiffusion.md` "Quality gates" section.

## RFdiffusion seed threading is verified per pod, not assumed

Trap #14. `workflow/scripts/run_rfdiffusion.py` ships with
`_SEED_THREADING_MODE = "per_design"` as the safe default. RFdiffusion's
`inference.num_designs=N inference.random_seed=base` *may* seed designs
at `base, base+1, ..., base+N-1` -- but this depends on the version of
the cloned repo and is not documented as a stable contract. **Fix**: on
each fresh pod or after a `git pull` in `/workspace/RFdiffusion/`, run
the five-minute recon at the top of the `run_rfdiffusion.py` docstring;
flip the constant to `"single_subprocess"` only if the three-design
smoke test confirms threaded seeding. The per-design fallback always
works but pays ~30-60 min of model-load overhead per 200-design batch.

## Never commit absolute paths in fixture YAMLs

Trap #16. Symptom: `pytest -q` is green locally but CI fails with
`FileNotFoundError: /home/user/.../tests/fixtures/.../something.pdb`.
Root cause: a fixture generator (e.g. `tests/fixtures/stage1/_make_fixtures.py`)
wrote an absolute path into a committed YAML via `str(some_absolute_path)`.
The path lives only in the generator's sandbox; once the repo is cloned
elsewhere, the path no longer resolves.

**Fix**: at *generation* time, write filesystem paths in committed
fixtures as either a bare filename (relative to the manifest's own
directory) or a repo-relative path -- never absolute. At *consumer*
time, resolve relative paths defensively: try
`(manifest_path.parent / relative).exists()` first, then fall back to
treating the path as cwd-relative. See
`workflow.scripts.run_rfdiffusion._resolve_cleaned_pdb`.

**Recurrence guard**: `tests/test_no_absolute_paths_in_committed_fixtures.py`
walks `tests/fixtures/`, parses every YAML, and asserts no string value
starts with `/`. Sub-100ms; catches the next leak before push.

**Why locally-clean isn't enough**: the failure mode requires running
tests from a different absolute path than the generator. Always verify a
PR in a fresh-checkout dir before push (`git clone . /tmp/clean && cd
/tmp/clean && uv sync --extra dev && uv run pytest -q`).

## RFdiffusion checkpoint integrity via pinned sha256

Trap #15. `bootstrap.sh --with-rfdiffusion-weights` downloads
`Complex_base_ckpt.pt` over `http://files.ipd.uw.edu/...` (the URL the
RFdiffusion README publishes; an `RFDIFF_CKPT_URL` env var lets the
operator override to https on the pod if it's available). Integrity is
verified by `sha256sum -c` against a pinned `RFDIFF_CKPT_SHA256`. **The
pin is empty in the initial commit by design**: the first
`--with-rfdiffusion-weights` run prints the computed sha256 and aborts
asking the operator to pin it. Subsequent runs verify. Same flow applies
in cycle 3 if the upstream checkpoint rotates -- unset, re-run, capture,
re-pin. No temporary-code ritual.
