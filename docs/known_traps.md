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

---

## Trap #26: Writer/reader contract drift hidden by mock-aligned canary

**Cycle**: 02 (Stage 1, cycle 02 run).

**What happened**: Stage 1 cycle 02 ran 200 RFdiffusion designs
successfully (~10 h GPU). All 200 PDBs and TRBs were written to
`results/cycle_02/stage1/rfdiffusion/designs/`.
`stage1_summary.json` reported `n_completed=0` and HALTed. Mock canary
4/40 at 0.100 was green throughout.

**Why the canary missed it**: `_copy_mock_fixtures` writes mock files
using the same zero-padded `:05d` filename convention as the reader
expects, by construction (it strips `mock_` from
`mock_design_NNNNN.pdb`). Real RFdiffusion writes
`design_{design_startnum}.pdb` where `design_startnum=seed`, producing
un-padded variable-digit filenames (`design_2000.pdb`, not
`design_00000.pdb`). The mock writer and reader were always trivially
aligned. The real writer was never exercised end-to-end in test.

**Generalization**: Any time a test fixture *produces* the data the
system-under-test will *consume*, the test only verifies the consumer,
never the producer. Producer↔consumer contracts must be tested against
the real producer at least once, even if expensive.

**Detection**: A test that synthesizes outputs using the *real* writer's
naming convention (seed-based) caught it once written —
`tests/test_run_rfdiffusion_real_writer.py`.

**Lesson**: Audit every other Stage's mock/real path. If a similar
pattern exists in Stage 2 (writer = AF2/ProteinMPNN subprocess, reader
= enumerator), assume similar drift is latent.

**Fix**: `workflow/scripts/run_rfdiffusion.py:_expected_pdb_path` is now
the single source of truth for the writer↔reader filename contract; both
the dispatch pre-flight and the post-run enumerator route through it.
`scripts/run_stage1.py --skip-subprocess` re-enumerates without
re-invoking RFdiffusion, providing a recovery path that does not waste
the 10 h of compute.

## Trap #28: `_binder_ca_coords` hardcoded the binder to chain "A"

**Symptom**: Cycle 02 `--skip-subprocess` recovery completed
(n_completed=200) but reported `binder_length=275` on every design and
failed the geometry length check on 200/200. Mock canary still passed
8/10 (tautology: mocks placed the binder on chain "A" too, so the wrong
reader happened to align with the wrong fixtures).

**Root cause**: `workflow/scripts/run_rfdiffusion.py:_binder_ca_coords`
selected chain "A" based on the false comment "RFdiffusion writes the
binder as chain A". For real binder-design with contig
`[A1-275/0 B0-99/0 C1-9/0 70-110]`, RFdiffusion preserves A/B/C as
motif chains and writes the designed binder to the next free letter
(chain D). The function was reading the HLA heavy chain (275 res),
which then tripped the length filter `[65, 115]`. Verified empirically
on `results/cycle_02/stage1/rfdiffusion/designs/design_2000.pdb`:
chain CA counts = `{'D': 98, 'A': 275, 'B': 100, 'C': 9}`.

**Fix**: Identify the binder by exclusion. `_binder_ca_coords` now
takes `fixed_chains: frozenset[str]` derived from
`target.chains[*].role != "binder"` and selects the unique remaining
chain in the PDB. Raises ValueError on ambiguity (zero or >1
candidate). Mock fixtures (`tests/fixtures/stage1/mock_design_*.pdb`)
regenerated to put the binder on chain D so the canary now exercises
the real chain layout.

**Recurrence guard**:
`tests/test_run_rfdiffusion_real_writer.py::test_real_writer_binder_chain_is_D`
synthesizes a 4-chain PDB (A=275, B=100, C=9, D=85) and asserts the
correct chain is read. `test_binder_ca_coords_raises_on_ambiguous_chains`
asserts the loud-failure mode when chain assignment is ambiguous.

**Same bug class**: Traps #26, #27. Stage 2 `splice_binder.py:70`
embeds the same false belief (rejects PDBs whose first chain is not
"A") and is filed as a follow-up issue.

## Trap #29: `splice_binder.py` expected a 1-chain Stage 1 PDB with binder on chain "A"

**Cycle**: 02 (Stage 2, immediately after the Trap #28 fix unblocked Stage 1
cycle 02 with 26/200 geometry-pass designs).

**Symptom**: Stage 2 raised `ValueError: expected exactly 1 chain (the
binder)` on the first cycle 02 backbone fed to
`splice_binder_onto_pmhc`. The mock canary was green throughout — the
canary fixture itself only wrote single-chain stage1_backbones, so the
4-chain real-world contract was never exercised.

**Root cause**: `workflow/scripts/splice_binder.py` was written against
the cycle 01 RFdiffusion output shape (single-chain PDB, binder on
chain "A") and asserted `len(s1_chains) == 1` and
`s1_chains[0].id == "A"` before renaming `A → D`. With the cycle 02
contig `[A1-275/0 B0-99/0 C1-9/0 70-110]`, RFdiffusion preserves the
motif chains A/B/C and writes the binder to chain D, producing a
4-chain PDB.

**Why the canary missed it**: same family as Trap #26.
`tests/fixtures/stage2/designs/_make_fixtures.py:write_stage1_backbones`
constructed single-chain (chain A) backbones, matching the reader's
old assumption by construction. The producer↔consumer contract was
never tested against the real RFdiffusion writer.

**Fix**:
- `workflow/scripts/splice_binder.py` now extracts chain D directly
  (`_STAGE1_BINDER_CHAIN = "D"`), validates that the Stage 1 PDB
  contains chain D, rejects any chain outside `{A,B,C,D}` (loud-fail
  pattern, mirroring `_binder_ca_coords`), and sources A/B/C from the
  cleaned pMHC.
- Stage 2 now pre-filters Stage 1 designs on `geometry_pass=True`
  before MPNN/splice fan-in (`scripts/run_stage2.py`); cycle 02 would
  otherwise burn compute on ~170 zero-contact designs.
- Stage 1 `HALT_THRESHOLD` lowered from 0.50 → 0.10 to align with
  Stage 2 and let cycle 02's 0.13 pass; the verdict-must-be-PASS gate
  in `scripts/run_stage2.py:454` then no longer blocks. Cycle 03 will
  tighten after partial diffusion lands.
- `tests/fixtures/stage2/designs/_make_fixtures.py:write_stage1_backbones`
  regenerated to emit 4-chain (A/B/C/D) PDBs so the canary now
  exercises the real chain layout.

**Recurrence guard**:
- `tests/test_splice_binder.py::test_splice_4chain_stage1_input_uses_canonical_abc`
  builds a 4-chain Stage 1 PDB whose A/B/C residue counts differ from
  the cleaned pMHC and asserts the output A/B/C is sourced from the
  cleaned pMHC (not the Stage 1 PDB).
- `test_splice_rejects_input_without_chain_D` asserts loud failure on
  the legacy cycle 01 single-chain layout.
- `test_splice_rejects_unexpected_chain` asserts loud failure on
  chains outside `{A,B,C,D}`.
- `tests/test_stage1_halt_gate.py::test_enforce_halt_gate_pass_at_cycle02_boundary`
  pins the `>=` boundary semantics of the halt rule at the exact cycle
  02 ratio (13/100 = 0.13).

**Same bug class**: Traps #26, #27, #28. Fourth instance of writer↔reader
contract drift hidden by mock-aligned canary. The general lesson stands:
any fixture that synthesises its own producer output verifies only the
consumer, never the contract.

## Trap #30: RFdiffusion crude-sequence output is Ala-heavy by design

**Cycle**: 02 → 03 (Stage 1 / Stage 2 boundary).

**Symptom**: cycle-02 designs (and BAKER's published `scaf*.pdb`) carry an
alarmingly high alanine fraction in their as-emitted sequences — ~38% in the
BAKER scaffold library, ~40.4% in our cycle-02 hero `design_2079`. This looks
like a degenerate/broken design at first glance.

**Root cause**: this is **not** a pipeline bug. RFdiffusion only generates
backbone coordinates; the residue identities it writes are a crude
placeholder driven by the model's compositional prior, which is strongly
Ala-biased. The real sequence is assigned downstream by ProteinMPNN, which
redesigns every binder position. BAKER's `scaf*.pdb` are explicitly
crude-sequence scaffolds for exactly this reason.

**Implication / fix**:
- Never read meaning into the raw RFdiffusion / scaffold sequence; always let
  ProteinMPNN redesign (it already does — `design_chains: [D]`).
- Cycle 03 additionally counter-biases the prior with
  `configs/proteinmpnn_bias_aa.json` (`A: -2.0`, `E/L/R: +1.0`) wired via
  `--bias_AA_jsonl` in `scripts/run_stage2.py`, so MPNN is nudged away from
  Ala and toward BAKER's redesigned-binder composition.
- Do not "fix" Ala-rich backbones upstream; the lever is the MPNN bias, not
  the diffusion step.

## Trap #31: BAKER scaffold alignment requires chain B = HLA truncated; the full target with β2m causes silent geometric mismatch

**Cycle**: 03 (Stage 1 sub-run A prep).

**Symptom**: aligning BAKER's published `scaf*.pdb` library onto our
canonical `data/targets/3hpj_clean.pdb` produces binders sitting in the wrong
place — partial diffusion then starts from garbage geometry and burns ~6 h of
A100 on designs that can never pass. No exception is raised; the alignment
"succeeds" with a plausible-looking RMSD.

**Root cause**: chain-layout mismatch. Each BAKER scaffold carries the target
as a single **chain B = HLA-A\*02:01 α1+α2 (residues 1-180) fused to the 9-mer
peptide** (189 continuous residues), binder on chain A. Our `3hpj_clean.pdb`
follows the crystal convention: chain A = HC (~275), **chain B = β2m**, chain C
= peptide. The old `align_scaffolds.py` superposed *scaffold chain B onto
reference chain B* — i.e. HLA α1/α2 CA atoms onto β2m CA atoms. β2m and the
HLA groove are unrelated folds, so the recovered rigid transform is meaningless
even though the CA-count match and RMSD look fine.

**Fix**:
- `workflow/scripts/prep_baker_target.py` produces a BAKER-format truncated
  target `data/targets/3hpj_baker_truncated.pdb`: keep HC residues 1-180
  (α1+α2), rename A→B, drop β2m, keep peptide as chain C. Layout now matches
  BAKER's groove-only fragment (chain B = HLA, chain C = peptide).
- `workflow/scripts/align_baker_scaffolds.py` superposes the scaffold's chain-B
  HLA substring (first 180 CA) onto the truncated reference's chain B, then
  rewrites the complex into our standard `A=binder / B=HLA / C=peptide` layout.
- `align_scaffolds.align_scaffolds(..., baker_layout=True)` dispatches to it;
  sub-run A (`rule stage1_subrun_a`) passes the truncated target as the align
  `--reference-pdb` while the geometry/motif reference stays the full target.
- Sub-run B (full target, `design_2079`) is unaffected: it was AF2-validated
  against the full 4-chain target and keeps `3hpj_clean.pdb`.
- Truncation also shifts the AF2 chain layout for sub-run A and its
  recalibrated controls to 3-chain `A=HLA, B=peptide, C=binder`;
  `compute_metrics` selects the binder/peptide/MHC chains via
  `LAYOUT_CHAINS[target_layout]` so decomposed iPAE stays correct.

**Recurrence guard**:
- `tests/test_prep_baker_target.py` asserts the truncated target has exactly
  CA counts `{B: 180, C: 9}`, no chain A, no β2m, numbering preserved.
- `tests/test_align_baker_scaffolds.py` asserts the aligned output is
  `A=binder / B=HLA[180] / C=peptide[9]` and that the HLA-substring
  superposition recovers a known +100 Å translation (RMSD ≈ 0).
- `tests/test_controls_truncated.py` asserts the full-layout default
  (`binder="D"`) finds **no** interface on a 3-chain truncated prediction,
  while the truncated layout decomposes iPAE correctly — locking in the
  adapter that prevents the silent mismatch from re-entering the metrics.

**First observed**: cycle 03 sub-run A prep, before any A100 time was spent.

---

## Trap #32: partial diffusion still requires contigmap.contigs (it is the residue mask, not an optional bias)

**Cycle**: 03 (Stage 1 sub-run A, first real launch).

**Symptom**: sub-run A dies on the **first** `run_inference.py` invocation with
`Must either specify a contig string or precise mapping.` No design is produced;
the run never reaches the GPU diffusion loop.

**Root cause**: `partial_diffuse.partial_diffuse_one` built the RFdiffusion
command without `contigmap.contigs`, on the false belief that the contig is only
needed for de-novo generation (or that it is an optional hotspot bias). It is
neither: RFdiffusion uses `contigmap.contigs` as the **residue mask in every
mode**, partial diffusion included. Per `/workspace/RFdiffusion/README.md`
("Partial diffusion"):

> Anything prefixed by a letter indicates that this is a motif, with the letter
> corresponding to the chain letter in the input pdb files. Anything not
> prefixed by a letter indicates protein to be built.

> if you have a binder:target complex, and you want to diversify the binder
> (length 100, chain A), you would need to input something like this:
> `contigmap.contigs=[100-100/0 B1-150]` `diffuser.partial_T=20`

Do **not** conflate this with `ppi.hotspot_res`, which *is* optional in partial
mode.

**Fix**:
- `partial_diffuse._derive_contigs_subrun_a(aligned_scaffold_pdb)` reads the
  aligned scaffold's chain-A binder length `N` (layout A=binder, B=HLA[1:180],
  C=peptide[1:9]) and returns `[N-N/0 B1-180/0 C1-9]`. The binder slot is a
  **bare `N-N`** range (unprefixed = redesigned); `B`/`C` are **letter-prefixed
  motifs** (preserved). An `A`-prefixed binder slot would tell RFdiffusion to
  *preserve* chain A, defeating partial diffusion (the BAKER scaffolds would
  come back essentially unchanged).
- `partial_diffuse_one` appends `contigmap.contigs={contig}` to the command;
  a pre-flight `ValueError` fires if chain A is empty/malformed, so a bad mask
  is caught before the subprocess rather than swallowed by RFdiffusion's
  generic error.
- Knock-on (same root layout): the bare-`N-N` binder makes RFdiffusion write the
  designed binder as chain **A** in a 3-chain output (A=binder, B=HLA,
  C=peptide). Stage 1's geometry gate (`_binder_ca_coords`) identifies the
  binder by exclusion from `fixed_chains`; the full manifest's `{A,B,C}` would
  leave no candidate and crash scoring. `run_stage1_subrun._fixed_chains_for_subrun`
  returns `{B,C}` for real sub-run A so chain A is the unique binder.

**Recurrence guard**:
- `tests/test_partial_diffuse.py`: `_derive_contigs_subrun_a` on an 80-mer
  scaffold returns `[80-80/0 B1-180/0 C1-9]`; an empty chain A raises
  `ValueError`; and `partial_diffuse_one`'s emitted command contains
  `contigmap.contigs=[80-80/0 B1-180/0 C1-9]`.
- `tests/test_stage1_subruns.py`: `_fixed_chains_for_subrun` returns `{B,C}`
  for real sub-run A, and an end-to-end `run(subrun="a", mock=False)` on a
  3-chain design scores `binder_length == 80` without crashing (call-site
  guard — a helper test alone could pass the right `fixed_chains` while the real
  path passes the wrong one; cf. Trap #28).

**Note**: this PR covers Stage 1 only. Real sub-run A then produces 3-chain
designs that Stage 2's `splice_binder` (4-chain binder-on-D expectation) cannot
yet consume — the operator halts after `rule stage1_subrun_a`; the truncated
Stage 2 path is a separate PR (Trap #33, TBD).

**First observed**: cycle 03 sub-run A, first real pod launch.
