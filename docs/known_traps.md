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

## Stage 1 binder lives on chain "A" -- always rename to "D" before AF2

Trap #17. Symptom: AF2 returns an iPAE of garbage (often <2 Å or >40 Å)
for spliced complexes despite ColabFold reporting a high ranking
confidence. Compute_metrics happily produces a single number; the halt
gate may even pass on individual cases. Root cause: Stage 1 writes each
RFdiffusion-generated binder as a single-chain PDB with chain id "A".
If that PDB is fed directly to ColabFold alongside an A/B/C cleaned
pMHC, the multimer assembly gets two chain "A"s collapsed into one
(ColabFold uses chain ids to determine the asym table), and the binder
is silently treated as part of the heavy chain. The iPAE then becomes a
within-chain distance metric -- meaningless.

**Fix**: `workflow/scripts/splice_binder.py` always renames the Stage 1
chain "A" to "D" while composing the 4-chain complex (A=HC, B=beta2m,
C=peptide, D=binder) and writes per-chain residue numbering starting at
1. The Stage 2 designs orchestrator (`scripts/run_stage2.py`) routes
every Stage 1 PDB through the splice helper -- never pass them directly
to ColabFold.

**Recurrence guard**: `tests/test_splice_binder.py::test_splice_chain_renaming`
+ `test_splice_rejects_input_with_chain_D` assert that the output
contains exactly A/B/C/D and that the splice rejects malformed inputs.

## Trap #18: bootstrap.sh had the wrong URL hash for `Complex_base_ckpt.pt`

**Symptom**: `bash bootstrap.sh --with-rfdiffusion-weights` reached step 5/7, attempted `wget http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Complex_base_ckpt.pt`, returned exit code 8 (server error), left a 0-byte file at the destination. Re-running just re-attempts the same broken URL.

**Why**: The hash segment `6f5902ac237024bdd0c176cb93063dc4` belongs to `Base_ckpt.pt`, not `Complex_base_ckpt.pt`. The IPD server returns 404 for the mismatched path. wget's default behavior on 404 is silent unless run with `--show-progress` or `-v`. The correct hash for Complex_base_ckpt.pt is `e29311f6f1bf1af907f9ef9f44b8328b` — verified against RFdiffusion's official `scripts/download_models.sh`.

**Fix**: Hotfix commit `fix(bootstrap): correct Complex_base_ckpt.pt URL hash + pin sha256`. URL corrected in `bootstrap.sh`. Sha256 pinned to `76e4e260aefee3b582bd76b77ab95d2592e64f00c51bf344968ab9239f3250bc`. Future pods don't hit this.

**Diagnostic shortcut**: Any time `bootstrap.sh`'s weight-download step exits with code 8 AND leaves a 0-byte file, the URL is wrong (not a network glitch). Verify with `curl -I <url>` before retrying. A 404 response means fix the URL, not the network.

**First observed**: cross-region pod migration, 2026-05-18.

---

## Trap #19: Pod volume vs Network volume confusion in RunPod's new deploy UI

**Symptom**: A pod was deployed with only a Pod volume (e.g., 50 GB) mounted at `/workspace`, no existing network volume attached. The Pod details panel shows "Pod volume, <size>" instead of "Network volume, <size>". Bootstrap would work, but (a) data is lost on pod termination, and (b) the Pod volume is typically too small for the full pipeline (RFdiffusion + ColabFold weights + 200 designs + AF2 outputs need ~30 GB minimum, with no headroom).

**Why**: RunPod's new deployment UI auto-creates a Pod volume by default if no Network volume is selected. The deployment panel's right-side summary toggles between "Network volume will be created" (= new empty volume, bad if you have an existing one) and "Network volume will be attached" (= the existing one). Easy to miss when configuring quickly.

**Verification protocol for future pod deployments**:
1. In the new RunPod deploy UI, locate the **Network volume** dropdown in the Compute section (between search box and Filter).
2. Click it → select the existing network volume by name.
3. Verify the right-side panel changes to show the volume name and "will be attached".
4. After deployment, immediately verify in the Details tab: should show "Network volume" with the correct size and name. If it shows "Pod volume" — wrong, terminate and redeploy.

**First observed**: cross-region pod migration, 2026-05-18 (pod `attractive_olive_moth`, AP-IN-1, terminated and replaced).

---

## Trap #20: Hardware drift in AF2 numerical outputs across GPU architectures is bounded but real

**Symptom**: Running the same pipeline (same code, same weights, same sha256-pinned dependencies, same inputs) on different GPU architectures produces non-identical AF2-multimer outputs. iPAE drifts by ~±0.5 Å, ipLDDT by ~±5 points, BSA by ~±100 Å² between A100 and H100. The rank-001 model/seed winner can differ.

**Why**: AF2-multimer ensembles 5 models per prediction; the rank-001 winner is decided by predicted-confidence margins tight enough to be flipped by floating-point non-determinism between SM86 (A100) and SM90 (H100). cuDNN heuristics and TF32/BF16 paths also differ. The differences are bounded but real — not a bug, intrinsic to the AF2 ensemble selection process.

**Implication for halt gates**:
- Bake ≥ 1 Å margin into iPAE-based gates.
- Bake ≥ 5 ipLDDT points into ipLDDT-based gates.
- For cross-pan ΔiPAE in Stage 3, the intrinsic noise floor from model selection alone is ~0.3–0.5 Å — design discrimination thresholds need to be > that.

**Validation evidence**: cycle 01 controls re-run, 2026-05-18, A100 SXM US-CA-2 → H100 SXM US-NE-1. Pos↔Neg iPAE gap held at 20.0 Å (threshold ≥10); ipLDDT gap at 59.8 (threshold ≥30). All five controls passed on both hardware. Worst-case positive drift: P2 iPAE +0.34 Å, N2 ipLDDT −4.4 points. See `results/cycle_01/stage2/metrics.json` (H100, canonical) vs `results/cycle_01/stage2/metrics_A100_baseline.json` (A100, sidecar).

**Practical rule**: hardware drift in AF2 numerical outputs is expected and bounded. Cycle 1 had ~2 Å of margin on every halt-relevant dimension; this is the design target for future thresholds.

---

## Trap #21: Skip `--use-gpu-relax` during Stage 2 AF2 triage

**Symptom**: `colabfold_batch` spends 30–50 % of its wall time on AMBER
relaxation that does not change iPAE/ipLDDT ranking on binder triage runs.
At cycle 02's 50-prediction fan-in this can add ~45 min of pure overhead
without affecting which designs survive the halt gate.

**Why**: Baker lab's `pmhc_fold.py` runs `do_relax=False` during
high-throughput screening for exactly this reason — the relax step is a
post-hoc geometry polish, useful for final structure inspection but
noise vs. ranking signal during triage.

**Fix**: `workflow/scripts/run_colabfold.py:run_real` and the inline
`_run_real_colabfold` in `scripts/run_stage2.py` both omit
`--use-gpu-relax`. Apply AMBER relax only to the cycle 03+ top survivors
if BSA precision becomes a bottleneck.

---

## Trap #22: ProteinMPNN `.fa` files start with the original input sequence

**Symptom**: Stage 2 aggregator emits `n_seqs_per_backbone + 1` records
per design instead of `n_seqs_per_backbone`, contaminating the fan-in
ranking with a non-designed sequence at the top.

**Why**: Upstream `dauparas/ProteinMPNN` writes each `.fa` file with the
fixed-scaffold input as the first record (e.g. header containing
`fixed_chains=['A','B','C']`, no `sample=` field), followed by the
sampled designs (`T=0.1, sample=1, ...`, etc.). The first record is for
provenance, not a design.

**Fix**: `workflow/scripts/aggregate_mpnn_outputs.aggregate` always
skips `fasta_records[0]` and starts at `[1:]`. The unit test
`test_aggregate_skips_scaffold_record` is the recurrence guard.

---

## Trap #23: ProteinMPNN outputs the **full complex** sequence per record

**Symptom**: Stage 2 sequences.jsonl carries 200+-residue sequences for
each "design" when binder length is supposed to be 70–110.

**Why**: ProteinMPNN's `.fa` body contains the entire input PDB's
chain-concatenated sequence in PDB chain order, regardless of which
chains were fixed vs. designed. Chains A/B/C are returned verbatim
(native sequence); only chain D contains designed residues. The
aggregator must slice the binder portion.

**Fix**: `workflow/scripts/aggregate_mpnn_outputs.aggregate` accepts
optional `binder_length` per design (joined from Stage 1's
`designs.jsonl`) and slices `seq[-binder_length:]`. The unit test
`test_aggregate_slices_binder_when_binder_length_provided` is the
recurrence guard.

---

## Trap #24: ColabFold multimer FASTA is colon-separated within a single record

**Symptom**: `colabfold_batch` silently treats a multi-record FASTA as
independent monomer inputs; the resulting per-id PDBs have only one
chain and iPAE is undefined / garbage.

**Why**: ColabFold's multimer convention differs from DeepMind's original
AF2-Multimer format. ColabFold expects exactly one `>id\n` header per
record, with chains in the body separated by `:`. No trailing colon.
Example: `>design_00042_seq01\nGSHSM...:IQRTPK...:RMFPNAPYL:MAEEL...\n`.

**Fix**: `scripts/run_stage2.py:_build_multimer_fasta` writes a single
record with `":".join((A, B, C, binder))` and asserts no trailing colon
in `tests/test_run_stage2.py::test_build_multimer_fasta_colon_separated_single_record`.

---

## Trap #25: `parse_multiple_chains.py` reads every file under `--input_path`

**Symptom**: ProteinMPNN errors mid-batch with "no atoms found" or
similar, or includes stale PDBs from prior cycles in the run.

**Why**: Upstream's helper script does a directory walk — no extension
filter, no manifest. Anything under `--input_path` (Snakemake markers,
`.ipynb_checkpoints/`, accidentally-copied unrelated PDBs) is parsed
and either crashes the script or pollutes the JSONL.

**Fix**: `scripts/run_stage2.py:_stage_input_pdbs` always symlinks the
intended spliced PDBs into a fresh `input_pdbs/` directory under the
per-cycle MPNN work dir, then points `--input_path` at it. Idempotent:
the directory is removed and recreated on each call.
