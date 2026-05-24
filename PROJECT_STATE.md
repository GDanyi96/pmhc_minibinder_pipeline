# pMHC-I Minibinder Pipeline — Project State

**This is the single authoritative state document.** Read this first in any new Claude thread.

*Last updated: Sunday, May 24, 2026 (post cycle-03 Stage 1 sub-run A complete: 72/150 geometry-pass = 48%, 3.7× cycle 02 baseline).*

---

## 1. Thirty-second status

- **Goal**: production-grade de novo pMHC-I minibinder design pipeline mirroring `HADRUP_JENKINS_2025` (RFdiffusion → ProteinMPNN → AF2 → MD → IMPAC-T cells). PhD application portfolio targeting Jenkins/Hadrup lab (DTU/Glasgow).
- **Repo**: `github.com/GDanyi96/pmhc_minibinder_pipeline`. Active branch: `feat/cycle-03-baker-truncation` (PR #15 DRAFT, contains contig fix + chain layout reconciliation, NOT merged to main).
- **Compute**: RunPod A100 SXM 80GB pod `925f65a88d93` on network volume in **US-WA-1**. Old US-CA-2 pod stopped (cold archive).
- **Workflow**: Claude Code Web for repo writes; RunPod web terminal for compute. Never push from pod.
- **Cycle status**:
  - Cycle 1 controls (full target): validated.
  - Cycle 2 Stage 1 + Stage 2: complete. Hero `design_2079_seq00` (99aa 4HB, iPAE 6.41, ipLDDT 91.07).
  - **Cycle 3 controls (truncated target): validated** — P1 iPAE 3.33, P2 3.73, N1 22.0, ~18 Å dynamic range. Halt thresholds recalibrated to iPAE≤6.0, ipLDDT≥92.
  - **Cycle 3 Stage 1 sub-run A (BAKER scaffold library + partial diffusion): complete, 72/150 (48%) geometry-pass, median 6.5 hotspot contacts, 98 min wall on A100.**
  - Cycle 3 Stage 1 sub-run B (hero seed partial diffusion): **deferred, blocked** — see §16.
  - Cycle 3 Stage 2 truncated path: **implemented (mock-validated), pending pod real-run** — see §16. Trap #33 filename fix landed in the same PR.
- **Immediate next step**: pod real-run of Stage 2 sub-run A (`snakemake results/cycle_03/stage2/subrun_a/stage2_summary.json --config cycle=03 mock=false -j1`) to produce final heroes. Then optionally sub-run B (needs `_derive_contigs_subrun_b`, out of scope of this PR).

---

## 2. Working compute

| Property | Value |
|---|---|
| Pod hostname | `925f65a88d93` |
| GPU | A100-SXM4-80GB |
| Region | US-WA-1 |
| Driver | 550.127.05 (max CUDA 12.4) |
| Network volume | mounted at `/workspace`, 200+ TB free |
| Container overlay disk | 30 GB (small, fills fast — see Trap #39) |
| Old pod | US-CA-2 `kk519bspkk` stopped, cold archive ~$7/mo |

### Two coexisting Python environments

| Env | Location | Manager | Purpose | Critical pins |
|---|---|---|---|---|
| Main pipeline | `/workspace/pipeline/.venv` | uv | ColabFold, ProteinMPNN, Snakemake, metrics, AF2 fan-in | JAX 0.4.34, torch 2.12+cu130, **nvidia-cudnn-cu12==9.1.0.70** (forced down from 9.22 to match driver 550) |
| RFdiffusion | `/workspace/miniconda3/envs/SE3nv` | conda + manual pip overlay | RFdiffusion inference only | Python 3.9, **torch 1.9.0+cu111**, dgl 0.9.1.post1, hydra-core, e3nn==0.3.5, sympy==1.11.1, mpmath, opt-einsum, opt_einsum_fx, dllogger (NVIDIA git), se3-transformer (from RFdiffusion/env/SE3Transformer), rfdiffusion (editable from /workspace/RFdiffusion). **Every pip install MUST use `--no-deps`** — transitive deps will pull torch 2.x and destroy env (Trap #36). |

`RFDIFFUSION_PYTHON=/workspace/miniconda3/envs/SE3nv/bin/python` set in pod shell. Orchestrator honors it.

### MANDATORY environment exports before any pod work

```bash
export TMPDIR=/workspace/.cache/tmp
mkdir -p $TMPDIR
```

Without this, RFdiffusion/pytest/ColabFold dump scratch into 30 GB overlay disk and fill it within minutes (Trap #39). Add to `~/.bashrc` doesn't help between session restarts — must `export` in the current shell before launching anything.

### External repos + models on pod

- `/workspace/RFdiffusion/` at SHA `2d0c003df46b9db41d119321f15403dec3716cd9`
- `/workspace/pMHCI_binder_design/`, `/workspace/dl_binder_design/`
- `/workspace/models/RFdiffusion/Complex_base_ckpt.pt` (sha256 `76e4e260aefee3b582bd76b77ab95d2592e64f00c51bf344968ab9239f3250bc`)
- `/workspace/models/RFdiffusion/Complex_Fold_base_ckpt.pt`

---

## 3. Source paper anchors

| Anchor | What it is |
|---|---|
| `HADRUP_JENKINS_2025_pMHC_minibinder_landmark_Science` | THE blueprint paper. RFdiffusion + ProteinMPNN + AF2 + MD + IMPAC-T. NY-ESO-1 + neoantigens. PDB 9NNF. |
| `BAKER_LAB_2025_pMHC_binder_specificity_Science` | Baker co-submission. 11 pMHCIs across 4 HLAs. wt1-5 = 85 aa de novo. mage-513 = privileged 102-aa scaffold (target-agnostic). Peptide-centric arcing for placement constraint. **In silico funnel stays on truncated target throughout — β2m only reintroduced for crystallography.** |
| `HOUSEHOLDER_GARCIA_2025_helical_TCR_mimic_NYESO1_Science` | Garcia co-submission. 4-helix bundle TCR mimic. 9.5 nM. Hamming proteome scan for in silico specificity. |
| `MARES_IOANNIDIS_2025_diffusion_pMHC_peptide_library_benchmark` | UC Berkeley/CZ Biohub. Sequence-based predictors (NetMHCpan, MHCflurry, HLApollo) blind to structurally valid novel peptides. Use as baselines, not oracles. |
| `BENTZEN_HADRUP_2019_TCR_fingerprinting_NatBiotech` | DNA-barcoded MHC multimer specificity. Foundation for experimental cross-reactivity readout. |
| `HADRUP_ESMO_2025_minibinder_clinical_perspective` | Clinical framing for cover letter / elevator pitch only. |

Always anchor claims to one of these. Mark structural prediction outputs as predictions, not facts.

---

## 4. Locked decisions

| Decision | Value | Source |
|---|---|---|
| Primary target | WT1 RMFPNAPYL / HLA-A\*02:01 (PDB 3HPJ) | Matches Baker wt1-5 |
| Cycle 1 calibration anchor (P2) | NY-ESO-1 SLLMWITQC / HLA-A\*02:01 (Jenkins NY1-B04, 145 aa) | HADRUP_JENKINS_2025 |
| **Cycle 03 target geometry** | **Truncated: HLA chain B = residues 1-180 (α1+α2 + start of α3), peptide chain C = 9 res. NO β2m, NO α3 below residue 180.** | BAKER_LAB_2025 Materials & Methods |
| **Cycle 03 in silico funnel convention** | Stay on truncated geometry through RFdiffusion → ProteinMPNN → AF2. β2m never reintroduced. | BAKER_LAB_2025 |
| Cycle 2 binder length | 70–110 aa | Baker de novo cluster centered on wt1-5 (85 aa) |
| Cycle 2 num_designs | 200 | Project scope |
| Hotspots cycle 02 | C1, C4, C6, C8 (peptide) + A65, A66, A150, A155 (MHC α1/α2) | `configs/target_wt1_a0201.yaml` |
| Hotspots cycle 03 | Same 8 residues (peptide C1/C4/C6/C8 + HLA A65/A66/A150/A155); resolved per-design vs aligned scaffold coords | `configs/target_wt1_a0201.yaml` |
| **RFdiffusion contig format** | `[N-N/0 B1-180/0 C1-9]` — **bare** `N-N` for binder (unprefixed = redesigned); letter-prefixed `B1-180`, `C1-9` = preserved motif | `/workspace/RFdiffusion/README.md` "Partial diffusion" — Trap #32 |
| **Chain layout cycle 02 (full target, de novo)** | RFdiffusion output: A=HC, B=β2m, C=peptide, D=binder (binder is the diffused chain, gets next free letter) | Empirically verified cycle 02, Trap #29 |
| **Chain layout cycle 03 sub-run A (truncated, partial diffusion from BAKER)** | RFdiffusion output: **A=binder, B=HLA[1:180], C=peptide** (input scaffold chain IDs preserved; BAKER scaffolds have binder on A) | Empirically verified cycle 03 |
| **Chain layout truncated controls / Stage 2 FASTA** | `HLA[1:180]:peptide:binder` → AF2 output A=HLA, B=peptide, C=binder | Matches BAKER convention; controls baseline |
| **`LAYOUT_CHAINS["truncated"]`** | `binder="C", peptide=("B",), mhc=("A",)` | `workflow/scripts/compute_metrics.py` |
| **Stage 2 FASTA for sub-run A designs (NOT YET IMPLEMENTED)** | Must reorder from design's A=binder/B=HLA/C=peptide → FASTA `HLA:peptide:binder` to match `LAYOUT_CHAINS["truncated"]` after AF2 | Pending Stage 2 truncated PR (§16, Trap #33-adjacent) |
| iPAE definition | Symmetric mean cross-chain PAE binder↔target | `workflow/scripts/compute_metrics.py` |
| AF2 recycles | 6 baseline, 3 for cycle 02 PoC | Matches Jenkins/Baker |
| Halt rule, cycle 1 controls | Separation-based: pos↔neg iPAE gap ≥ 10 Å, ipLDDT gap ≥ 30 | Cycle 1 reformulation |
| Halt rule, Stage 1 cycle 2/3 | `fraction_geometry_pass ≥ 0.10` | Cycle 02 calibration data |
| **Halt rule, Stage 2 cycle 3 (NEW, truncated)** | **`halt_cut_ipae_max: 6.0`, `halt_cut_iplddt_min: 92.0`** (tightened from cycle 02's 10.0/88.0) | `metrics_truncated_baseline.json`: P1 iPAE 3.33, N1 22.0, ~18 Å dynamic range |
| Seed convention | `seed = cycle * 1000 + offset` per-stage in `configs/seeds.yaml` with reserved ranges | `configs/seeds.yaml` |
| Install model | Native, no Docker | Operational |
| MD validation | Cycle 4+ (deferred) | Project plan |
| Cross-pan off-targets (Stage 3) | MART-1 ELAGIGILTV + HIV KLTPLCVTL | Inherited P3 + Baker |

---

## 5. Pipeline state by stage

| Stage | Status | Key files |
|---|---|---|
| Stage 0 (target prep, full) | ✓ Working, cycle 02 used | `workflow/scripts/prep_target.py`, `configs/target_wt1_a0201.yaml`, `data/targets/3hpj_clean.pdb` |
| **Stage 0 (target prep, truncated)** | ✓ **NEW cycle 03** — strips β2m (chain A) and α3 (residues 181+ of HC), preserves coords | `workflow/scripts/prep_baker_target.py`, `data/targets/3hpj_baker_truncated.pdb` (B=180 / C=9). **Coordinate frame is preserved from full target** (verified: A65 CA in full = B65 CA in truncated = (30.574, -3.158, 7.961)) |
| **Stage 1 sub-run A (BAKER scaffolds + partial diffusion)** | ✓ **NEW cycle 03 — 72/150 geometry-pass = 48%.** 152 BAKER scaffolds aligned to truncated target, 150 partial-diffused with partial_T=15. Median hotspot contacts on passing: 6.5 (cycle 02 was ~0). **Symlink hack applied for filename recovery — Trap #33.** | `workflow/scripts/align_baker_scaffolds.py`, `workflow/scripts/partial_diffuse.py`, `scripts/run_stage1_subrun.py`, `configs/rfdiffusion_subrun_a.yaml`, `data/scaffolds/baker_library/` (152 scaffolds symlinked) |
| Stage 1 sub-run B (hero-stitched partial diffusion) | 🔄 **Stubbed/deferred** — needs `_derive_contigs_subrun_b` helper (current code's `_derive_contigs_subrun_a` assumes BAKER 3-chain layout, would crash on 4-chain hero PDB) + Trap #33 fix | `data/seeds/design_2079_binder.pdb` (parked), `configs/rfdiffusion_subrun_b.yaml` |
| Stage 1 (cycle 02 de novo path) | ✓ Working, 200 designs ran, 13% geometry-pass | `scripts/run_stage1.py`, `workflow/scripts/run_rfdiffusion.py` |
| Stage 2 — controls (full target) | ✓ Validated cycle 1 | `scripts/run_controls.py --target=full` |
| **Stage 2 — controls (truncated)** | ✓ **NEW cycle 03** — validated, recalibrated halt gates | `scripts/run_controls.py --target=truncated`, `results/cycle_03/controls_truncated_baseline/metrics_truncated_baseline.json` |
| Stage 2 — designs (cycle 02 full target) | ✓ Complete, design_2079_seq00 hero | `scripts/run_stage2.py`, `workflow/scripts/splice_binder.py` (4-chain expected) |
| **Stage 2 — designs (cycle 03 truncated sub-run A)** | ✓ **NEW cycle 03 — implemented (mock-validated; pending pod real-run).** Forked `splice_binder_subrun_a` (3-chain → A=HLA, B=peptide, C=binder); `run_stage2 --target-layout truncated` threads layout through splice/MPNN(`--chain_list C`)/FASTA(`HLA:peptide:binder`)/metrics. Snakemake rule `stage2_subrun_a` (explicit target). | `scripts/run_stage2.py`, `workflow/scripts/splice_binder.py::splice_binder_subrun_a`, `workflow/rules/02c_stage2_subrun_a.smk`, `tests/test_stage2_subrun_a.py` |
| Stage 3 (in silico cross-pan) | Architected, not coded | spec TBD |
| Stage 4 (embedding diversity curation) | Architected, not coded | spec TBD |
| Stage 5 (MD validation) | Cycle 4+ | deferred |
| Stage 6 (active learning) | After experimental round | deferred |

---

## 6. Operational rules

Non-negotiable. Each rule exists because skipping it cost hours.

1. **Never push from the pod.** Pod has no GitHub auth.
2. **Replace `git pull` with `git fetch origin && git reset --hard origin/<branch>` on the pod.** Pipeline writes to tracked files.
3. **All pod jobs ≥5 min run in tmux.**
4. **`uv run python` for the pipeline; SE3nv's python only for RFdiffusion inference.** Don't cross.
5. **`os.environ.pop("LD_LIBRARY_PATH", None)` at the top of every CUDA-touching script.**
6. **Mock tests use `--cycle 99`. NEVER `--cycle 01/02/03`.**
7. **All committed fixtures use relative paths.**
8. **Fresh-clone pre-push gate before any final push** that touches code or fixtures.
9. **Dependency-import smoke test before declaring a pod "READY".**
10. **Mock CI green ≠ real-run validated.** Always real-run against published controls.
11. **When stuck on env issues for >30 min, step back** and ask whether the framing is wrong.
12. **Sanity-check compute math before locking scope.**
13. **Verification block after multi-file operations** with expected outputs labeled inline.
14. **PR cross-check before merging anything that touches gates, fixtures, or shared contracts.** Standard 8-section diff check.
15. **New CC session per PR** to avoid context bloat and decision-drift.
16. **🆕 `export TMPDIR=/workspace/.cache/tmp` in every new shell before any pipeline/pytest/RFdiffusion work.** The 30 GB container overlay fills within minutes otherwise (Trap #39).
17. **🆕 Read the canonical source (paper, README, repo) BEFORE improvising fixes.** Cycle 03 launch lost hours improvising contig syntax instead of reading `/workspace/RFdiffusion/README.md` upfront (Trap #32).
18. **🆕 Switch strategy after 2–3 failed attempts at version-juggling deps.** Cycle 03 SE3nv torch-1.9 reconstruction was 4+ hours of cascading conflicts; the answer was always `--no-deps` discipline (Trap #36).
19. **🆕 Snakemake positional target args must precede `--config`.** `--config` consumes variadic until next flag, so `snakemake --config mock=false cycle=03 TARGET.json` fails with `ValueError: Invalid config definition`. Correct: `snakemake TARGET.json --config mock=false cycle=03 -j1` (Trap #34 sibling).
20. **🆕 Verify writer/reader filename agreement on EVERY new tool call site, not just at the cycle boundary.** Cycle 02 = Trap #28 (zero-padded vs seed-named). Cycle 03 = Trap #33 (single vs double-suffix). Same root failure mode, different cycle.

---

## 7. Trap book

Numbered for stable cross-reference. Authoritative numbering — `docs/known_traps.md` matches.

### Cycle 1 era — environment + code stack

1–12: (unchanged — see prior PROJECT_STATE for Docker/LD_LIBRARY_PATH/ColabFold flat-outputs/AF2 PAE rename/rank-based halt/mock fixture paths/.gitignore re-allow/region-locked volumes/A100 scarcity/JAX 0.4.34 pin/MSA-server vs local-DB/pod-no-GitHub-auth).

### Cycle 2 setup era — pipeline logic

13. **Geometry-pass thresholds can't be set without empirical data.**
14. **RFdiffusion seed-threading mode**: per-design subprocess is the safe default.
15. **Self-bootstrapping sha256 pin pattern.**
16. **Absolute paths in committed fixtures** cause CI/local divergence.
17. **~~Chain renaming Stage 1 → Stage 2~~** (factually wrong as written — corrected by Trap #29).

### Cycle 2 hardware migration era

18. **`bootstrap.sh` had wrong URL hash for `Complex_base_ckpt.pt`** — silent 404, zero-byte file.
19. **RunPod's new deploy UI** auto-creates a Pod volume by default. Verify Network volume.
20. **Hardware drift in AF2 across GPU architectures is bounded but real.** iPAE ±0.5 Å, ipLDDT ±5 between SM86 (A100) and SM90 (H100).

### Cycle 2 RFdiffusion env era

21. **`bootstrap.sh` does not install RFdiffusion's Python deps.**
22. **SE3nv canonical install yields CPU-only torch.**
23. **H100 (SM90) is incompatible with canonical SE3nv (CUDA 11.1 / torch 1.9).** Use A100 for RFdiffusion.
24. **Pod volume cannot be resized on a running pod.**
25. **`bootstrap.sh` writes target manifest to wrong path.**
26. **`UV_HTTP_TIMEOUT` default (30s) is too short for giant CUDA wheels.** Set 600.
27. **`configs/target_2bnr_a0201.yaml` has stale length range from cycle 1.**

### Cycle 2 writer/reader contract drift class (PRs #11–#13)

28. **Stage 1 enumerator filename convention drift** (PR #11). Real writer uses seed-named files; reader assumed zero-padded loop index. Mock fixtures used reader's wrong convention. Fix: centralize filename construction + `--skip-subprocess` re-enumeration flag.
29. **Stage 1 binder chain misidentification** (PR #12). Code claimed "binder is chain A"; RFdiffusion actually writes binder as next free chain letter (D in 4-chain layout). Fix: identify binder by exclusion, raise on ambiguity.
30. **Stage 1 Snakemake per-design output sentinel was hardcoded.**
31. **Stage 2 `splice_binder.py` 4-chain expectation** (PR #13). Same chain-A belief as Trap #29 one stage downstream. Fix: extract chain D from 4-chain Stage 1 PDB.

### 🆕 Cycle 3 era — partial diffusion + truncated target

32. **Partial diffusion STILL requires `contigmap.contigs` — it is the residue mask, not an optional bias.** Cycle 03 launch dies on first `run_inference.py` with `Must either specify a contig string or precise mapping.` Earlier guidance that the key was "unused by partial diffusion" was wrong. Per `/workspace/RFdiffusion/README.md` "Partial diffusion": `contigmap.contigs=[100-100/0 B1-150]` style — **bare `N-N` for binder (unprefixed = redesigned), letter-prefixed for motif (preserved).** Mistakenly using `A1-N` for binder slot would tell RFdiffusion to preserve chain A as motif → defeats partial diffusion entirely. **Fix**: `_derive_contigs_subrun_a` in `partial_diffuse.py` derives `[N-N/0 B1-180/0 C1-9]` per-design from aligned scaffold's chain-A length. Pre-flight `ValueError` guard if chain A absent. Don't conflate with `ppi.hotspot_res` (IS optional in partial mode). **Recurrence guard**: `test_derive_contigs_subrun_a_80mer`, `test_partial_diffuse_one_passes_contigs`.

33. **✅ RESOLVED (Stage 2 truncated PR) — Partial diffusion filename mismatch.** `partial_diffuse_one` constructed `inference.output_prefix=designs/design_{seed}` AND `inference.design_startnum={seed}`. RFdiffusion appends `_{startnum}` to whatever prefix is passed → produced `design_{seed}_{seed}.pdb` (double-suffix). Enumeration in `run_stage1_subrun.py` and `run_stage2.py` looks for `design_{seed}.pdb` (single). Cycle 03 sub-run A symptom: 150 designs successfully written, 0 enumerated, `subrun_summary.json` reports `n_completed: 0`, no error thrown. **Recovery hack (cycle 03 session)**: symlinks `design_NNNN.pdb -> design_NNNN_NNNN.pdb`, then re-trigger Stage 1 (existence check follows symlinks, skips RFdiffusion subprocess, runs only geometry scoring). **Permanent fix (LANDED)**: `partial_diffuse.partial_diffuse_one` now pins the prefix base to `out_prefix.parent / "design"` (caller stem ignored) and lets `design_startnum` produce the suffix → `design_{seed}.pdb` ✓. The symlink hack is no longer needed for new runs. **Tests**: `test_partial_diffuse_one_prefix_has_no_seed` (writer) + `test_subrun_a_writer_reader_filename_roundtrip` (real command + enumeration glob together). This is **Trap #28's twin** — same writer/reader filename drift class, different cycle, same lesson.

34. **`snakemake --config mock=false` stores "false" as STRING; `bool("false")` is True.** Snakefile uses `MOCK: bool = bool(config["mock"])` which evaluates `bool("false") == True` — pipeline silently runs in mock mode when user requested real. **Fix (applied pod-locally + committed)**: `_to_bool()` helper handling string + bool cases. Affects every Snakemake config key that the user passes as string flag.

35. **`uv venv` hardlinks files from `~/.cache/uv`.** `pip install --force-reinstall` and `uv pip install` both fail to replace pinned cuDNN wheel because the cache still hardlinks the old version. **Fix**: `uv cache clean` + `rm -rf .venv` + `uv sync --all-extras` + final `pip install --force-reinstall --no-deps nvidia-cudnn-cu12==9.1.0.70`. Cycle 03 spent ~2h on this before figuring out the cache layer.

36. **`pip install` in SE3nv without `--no-deps` ALWAYS pulls torch 2.x.** Every transitive dep modern enough to be on PyPI requires torch >= 2.0, which overwrites SE3nv's torch 1.9 → RFdiffusion immediately broken with cuDNN/CUDA errors. **Rule**: SE3nv installs use `pip install --no-deps PACKAGE` for everything, then manually resolve missing deps one at a time using `--no-deps` for each. Cycle 03 SE3nv reconstruction took ~4h because this rule wasn't enforced consistently. Packages currently installed (all `--no-deps`): se3-transformer (from /workspace/RFdiffusion/env/SE3Transformer), rfdiffusion (editable from /workspace/RFdiffusion), opt_einsum, opt_einsum_fx, e3nn==0.3.5, sympy==1.11.1, mpmath, dllogger (from NVIDIA git), icecream, decorator, jedi, hydra-core, omegaconf, pyrsistent, pynvml.

37. **RunPod network volumes are REGION-LOCKED.** Files written on US-CA-2 volume invisible to US-WA-1 pods. Cycle 03 hero PDB (`design_2079_seq00*.pdb`) had to be re-uploaded after pod migration. **Mitigation**: keep critical artifacts as git-tracked fixtures OR re-upload as part of pod-restart procedure.

38. **Stale mock outputs in `results/cycle_NN/` cause Snakemake to skip real stages.** If you `--config mock=true cycle=03 -j1`, then later `--config mock=false cycle=03 -j1`, Snakemake sees existing outputs and reports "Nothing to do." **Rule**: always `rm -rf results/cycle_NN/stageX results/cycle_NN/stageX+1` before launching a real run from mock state.

39. **Container overlay disk (30 GB) is TOO SMALL for /tmp during RFdiffusion + pytest.** Symptom: `OSError: [Errno 28] No space left on device: '/tmp/...'`. Pytest `tmp_path` fixture creates `pytest-0` accumulating subdirs that can't extend past 10 numbered tries. RFdiffusion writes Hydra logs to `/tmp` by default. **Fix**: `export TMPDIR=/workspace/.cache/tmp` in EVERY new shell before any pipeline work. Network volume `/workspace` has 200+ TB.

40. **🆕 HLA-CA RMSD is NOT predictive of BAKER scaffold transfer quality.** Counter-intuitive cycle 03 finding. Cross-tab: low-RMSD scaffolds (<1Å, the BAKER A\*02:01-native population, n=56) had 38% pass; high-RMSD scaffolds (≥3Å, the cross-allele population, n=96) had 54% pass. Reasoning: well-aligned scaffolds inherit their ORIGINAL peptide-targeting bias (designed against MART-1/gp100/NY-ESO etc., not WT1); badly-aligned scaffolds get "twisted" by the alignment routine which incidentally repositions binders. **Methodological consequence**: cycle 04 should filter by binder-to-peptide CA proximity, not HLA-HLA structural similarity. **Recurrence guard**: any pre-filter on scaffold library by HLA-CA RMSD is the WRONG mental model — log this explicitly when CC proposes such filters.

41. **🆕 Chain identities downstream of RFdiffusion must match the upstream tool's ACTUAL output convention — per sub-run, not per stage.** A tool that assumes "binder is always chain X" computes garbage silently when a new sub-run uses a different layout. Cycle-02 designs are 4-chain (binder=D); sub-run A designs are 3-chain (binder=A from RFdiffusion → C after the truncated splice). A metrics call hardcoded to `binder="D"` on a 3-chain truncated prediction finds no chain D → empty interface → iPAE `+inf`, no error. Same family as Traps #29/#31. **Fix**: the truncated Stage 2 path forks `splice_binder_subrun_a` (A=HLA, B=peptide, C=binder) and threads `target_layout` through `run_stage2` so splice / MPNN `--chain_list C` / FASTA `HLA:peptide:binder` / metrics all agree on `LAYOUT_CHAINS["truncated"]`; the cycle-02 4-chain path is untouched. **Rule**: document a new sub-run's chain layout end-to-end before writing code (see §12.i/j). **Recurrence guard**: `tests/test_stage2_subrun_a.py` (splice reorder + 4-chain rejection, FASTA order, MPNN chain dispatch, truncated metrics decomposition, mock end-to-end finite iPAE) + `tests/test_controls_truncated.py`.

42. **🆕 `stage2_subrun_a` mock branch bypasses its declared input (temporary layout mismatch).** The rule's **mock** branch reads dedicated 3-chain truncated fixtures (`tests/fixtures/stage2/designs_truncated/`) via `params.mock_summary`/`mock_manifest` instead of consuming `input.subrun_summary`. **Why**: the `stage1_subrun_a` mock still synthesizes **4-chain** designs (`_synthesize_mock_design` writes A/B/C+D), which `splice_binder_subrun_a` (needs 3-chain A=binder/B=HLA/C=peptide) correctly rejects — so the truncated mock path uses purpose-built 3-chain fixtures. The declared input still gates DAG ordering, and the **real** branch reads `input.subrun_summary` correctly. **Action**: re-couple the mock branch to `input.subrun_summary` once the `stage1_subrun_a` mock emits truncated 3-chain designs. A `TODO(Trap #42)` comment marks the mock branch in `workflow/rules/02c_stage2_subrun_a.smk`. **Recurrence guard**: `tests/test_stage2_subrun_a.py::test_subrun_a_stage2_mock_end_to_end`.

---

## 8. Pending CC work (queue, prioritized for cycle 03 continuation)

### 🔴 Priority 1 — required to continue cycle 03

1. **Stage 2 truncated path PR.** Splice for sub-run A (3-chain input, NOT 4-chain like cycle 02), FASTA construction in `HLA[1:180]:peptide:binder` order to match `LAYOUT_CHAINS["truncated"]`, MPNN designed-chain logic for chain-A binder, AF2 metric routing with `target_layout=truncated`. Tests: synthetic sub-run A design PDB → FASTA chain order matches expectation; mock cycle=99 integration covering sub-run A full Stage 2 path. New traps from this PR likely: chain-layout mismatch handling (#33-adjacent).

2. **Trap #33 filename fix PR.** Patch `partial_diffuse.partial_diffuse_one` to construct `out_prefix` without seed embedded; rely on `design_startnum` for suffix. Writer/reader contract test. Removes the symlink hack from cycle 03 recovery path. Affects both sub-run A (already-burned, recovered via symlinks) and sub-run B (would break the same way if launched).

3. **Sub-run B contig derivation PR** (only after #1 + #2 land). Add `_derive_contigs_subrun_b` helper for 4-chain hero-stitched layout (A=HC, B=β2m, C=peptide, D=binder). Dispatch from `partial_diffuse_one` based on input PDB chain composition or explicit sub-run argument.

### 🟡 Priority 2 — operational improvements

4. **Update `PROJECT_STATE.md`** (this doc) once committed, and update `docs/known_traps.md` with Traps #32–#40 in canonical form.

5. **Reproducibility hardening**: explicit pin file for SE3nv overlay packages (`requirements_se3nv_overlay.txt`) generated from current pod state so the env can be reconstructed without 4h debugging.

6. **Setup_cycle03_inputs.py predictions glob mismatch** (minor, surfaced during pod migration when re-uploading cycle 02 hero).

### 🟢 Priority 3 — backlog from earlier cycles

7. Commit drafted contract docs (`docs/contracts/rfdiffusion.md`, `proteinmpnn.md`, `alphafold2.md`).
8. Stage 2 mock canary audit for tautology (would catch Trap #33-class issues).
9. Trap #17 prose correction in `docs/known_traps.md`.
10. Restore A100 baseline as canonical `metrics.json`.
11. Renumber legacy topic-headed traps in `docs/known_traps.md`.
12. Extend `bootstrap.sh` to install miniconda + SE3nv conda env + GPU torch overlay.
13. Update `configs/target_2bnr_a0201.yaml` length range.
14. Parameterize controls expected-iPAE (currently hardcoded to full-target expectations; cycle 03 truncated controls run separately).

---

## 9. Standard fresh-pod sequence

*(See `docs/pod_quickstart.md` once committed. Key additions for cycle 03+: `export TMPDIR=/workspace/.cache/tmp` BEFORE any work; verify SE3nv torch is 1.9.0+cu111 via `python -c "import torch; print(torch.__version__)"`; verify `/workspace` has space via `df -h /workspace`.)*

---

## 10. Communication preferences

- Direct and technical. No "great question!"
- Design decisions: (1) recommended choice, (2) rationale anchored to a source paper, (3) known failure modes.
- Distinguish in silico evidence from experimental evidence.
- Code blocks, tables, structured lists. Prose only for synthesis.
- **One step at a time.** Wait for paste-back, then next step.
- **Verification block** after any multi-file operation.
- **Acknowledge own errors directly.** No apologies — correct and move on.
- For CC commits: user says "push directly to main" OR "push to feature branch, I'll merge via PR." Without that, CC defaults to feature branch.
- **Push back on incorrect claims with evidence**, not capitulation.
- **🆕 Always include bash code blocks for commands, never prose descriptions.** Multiple cycle 03 sessions wasted by Claude describing what to run instead of pasting the exact block.
- **🆕 Stop chasing version compatibility past 2-3 attempts; switch strategy.** Cycle 03 SE3nv reconstruction taught this — when transitive deps keep clashing, the answer is `--no-deps` discipline, not "try one more version combination."

---

## 11. Tool-boundary contracts (drafted, awaiting commit)

*(Unchanged — see prior PROJECT_STATE §11 for full RFdiffusion / ProteinMPNN / AF2 contract drafts. Cycle 03 adds an implicit additional contract for partial diffusion which should be appended to `docs/contracts/rfdiffusion.md` when committed: contig format, motif preservation in input chain IDs, partial_T schedule semantics.)*

---

## 12. What's fragile (read before touching adjacent code)

*(Sections a–h unchanged — see prior PROJECT_STATE §12.)*

### 🆕 i) Sub-run A vs sub-run B contig derivation is currently NOT polymorphic
`partial_diffuse.partial_diffuse_one` unconditionally calls `_derive_contigs_subrun_a`. Per CC's PR scope note, "deriving the sub-run-A contig unconditionally here is acceptable for this scope." Sub-run B would crash on this. When sub-run B path is implemented, dispatch must be explicit (subrun argument or input-PDB chain composition check). Don't let sub-run B silently fall through to the sub-run A helper.

### 🆕 j) `LAYOUT_CHAINS["truncated"]` is correct for controls and Stage 2 FASTA-builder OUTPUT, but **NOT for raw sub-run A design PDBs**
- Controls produce FASTA in `HLA:peptide:binder` order → AF2 output A=HLA/B=peptide/C=binder → matches `LAYOUT_CHAINS["truncated"]`.
- Sub-run A RFdiffusion outputs: A=binder/B=HLA/C=peptide (different convention because RFdiffusion preserves input scaffold chain IDs).
- These DON'T MATCH. Stage 2 must reorder sub-run A's chains during FASTA construction to convert A=binder/B=HLA/C=peptide → FASTA `HLA:peptide:binder`.
- Without this reorder, AF2 outputs would have binder on chain A and `LAYOUT_CHAINS["truncated"]` (which expects binder=C) would compute iPAE between wrong chains → garbage metrics.

### 🆕 k) Coordinate frame between `3hpj_clean.pdb` (full) and `3hpj_baker_truncated.pdb` (truncated)
**Verified preserved** (identical A65/B65 coords). This is required for the hotspot CA distance check to work (hotspot_xyz from full target manifest, design coords from truncated frame). If `prep_baker_target.py` ever re-orients or re-centers the truncated target, hotspot contacts become noise. **Recurrence guard**: any future Stage 0 change to truncation routine must include a coordinate-preservation assertion test.

---

## 13. Engineering principles surfaced (cycles 02–03 cumulative)

*(Cycle 02 principles a–g unchanged.)*

### 🆕 h) Cheap diagnostics buy expensive certainty
Cycle 03 launch session: a 30-second single-design RFdiffusion smoke caught the missing `contigmap.contigs` arg before the full 150-design run was launched. Single-design smoke is ~$0.05; full launch is ~$3-5. Always run a single-design smoke before any production GPU launch when the wrapper or contigs have changed.

### 🆕 i) Recovery affordances pay back per-instance, not just per-cycle
Trap #28's `--skip-subprocess` flag from cycle 02 PR #11 had a direct analog in cycle 03: the `if not out_pdb.exists()` natural existence check in `run_stage1_subrun.py` combined with symlinks let cycle 03's recovery use the same pattern. Each affordance built once is reused multiple times across cycles. Future stages should bake re-enumeration paths in by default.

### 🆕 j) "n_completed: 0" with healthy wall time is a writer/reader bug, not a design failure
Cycle 03 sub-run A initial summary reported `n_attempted: 150, n_completed: 0, wall_minutes_total: 98.5`. The instinct to read this as "designs failed" is wrong — 98 min of wall time with no errors thrown means designs WERE produced. n_completed=0 = enumeration couldn't find them = filename mismatch. **Diagnostic rule**: when summary counts say zero but wall time is reasonable, FIRST check `ls outputs_dir/` before assuming runtime failure.

### 🆕 k) Halt threshold migration must follow controls migration
Cycle 02 halt gate (iPAE≤10.0, ipLDDT≥88.0) was tuned against full-target controls. Cycle 03 truncated controls produced ~3 Å tighter iPAE distribution (P1: 4.5 full → 3.33 truncated; P2: 6.5 → 3.73). Migrating designs to truncated geometry WITHOUT migrating halt thresholds would make the halt gate too loose. **Rule**: whenever target geometry changes (truncated vs full, different HLA allele, etc.), recalibrate halt thresholds against the matching control panel BEFORE running designs.

### 🆕 l) Mental models can survive past their empirical anchor — confirm before relying
The "HLA-CA RMSD correlates with scaffold transfer quality" intuition felt obviously true (clean alignment = clean transfer). Cycle 03 data overturned it (Trap #40). Any intuition that hasn't been empirically confirmed in the current cycle's data should be flagged "untested" not "obvious" — especially when it would gate which scaffolds are pre-filtered out of expensive runs.

---

## 14. After cycle 03 sub-run A — queued work

1. **Stage 2 truncated path PR** (CC, Priority 1 item #1 above). Plumb the 72 designs through ProteinMPNN → AF2-multimer. Halt gate iPAE≤6.0, ipLDDT≥92.0. Targets to beat: P1 iPAE 3.33, P2 3.73; N1 22.0.
2. **Trap #33 filename fix PR** (CC, Priority 1 item #2). Removes symlink hack from recovery path.
3. **Optional Sub-run B** (CC, Priority 1 item #3, then GPU run). Hero seed partial diffusion. ~30 designs target.
4. **Cycle 04 plan** (artifact): scaffold pre-filter by binder-to-peptide CA proximity (NOT HLA-CA RMSD per Trap #40); length pre-filter to 70-110 (drops the 15 length-out-of-range failures up front); hotspot ablation experiment (with vs without `ppi.hotspot_res` for partial diffusion, paired on well-placed scaffolds); promote `mage-513` as privileged donor scaffold for a third sub-run.
5. **Cross-pan (Stage 3)** after Stage 2 truncated validates: MART-1 ELAGIGILTV + HIV KLTPLCVTL on the truncated convention.
6. **MD validation (Stage 5)**: cycle 04+.
7. **Application portfolio**: notebook draft anchored to:
   - Cycle 02: de novo baseline, 13% pass, hero design_2079
   - Cycle 03: scaffold partial diffusion 3.7× improvement (48%), HLA-RMSD non-predictiveness as methodological finding
   - Trap book as engineering-rigor evidence (Traps #28–#40 cumulative)

---

## 15. Cycle 02 numerical state (snapshot)

| Metric | Value |
|---|---|
| Cycle 02 Stage 1 designs attempted | 200 |
| Cycle 02 Stage 1 designs completed | 200 |
| Cycle 02 motif RMSD distribution | ~0.12–0.13 Å |
| Cycle 02 binder length distribution | min 70, median 89, max 110 |
| Cycle 02 Stage 1 geometry pass rate | 26/200 = 13% |
| Cycle 02 failure mode breakdown | 173/174 insufficient_hotspot_contacts, 4/174 internal_ca_clash, 3/174 multiple |
| Cycle 02 placement quality | 131/200 (66%) with ZERO hotspot contacts |
| Cycle 02 Stage 1 wall time | ~600 min (~10h) on A100 |
| Cycle 02 Stage 2 hero | design_2079_seq00 (99aa 4HB, iPAE 6.41, ipLDDT 91.07, 40.4% Ala) |
| Cycle 02 bugs caught + fixed during the run | 4 (Traps #28–#31), 4 PRs (#11–#13) |
| RFdiffusion git SHA used | `2d0c003df46b9db41d119321f15403dec3716cd9` |
| RFdiffusion Complex_base_ckpt.pt SHA256 | `76e4e260aefee3b582bd76b77ab95d2592e64f00c51bf344968ab9239f3250bc` |

---

## 16. 🆕 Cycle 03 numerical state (snapshot for portfolio)

### Truncated controls baseline

| Control | iPAE (Å) | ipLDDT | Verdict |
|---|---|---|---|
| P1 (Baker wt1, on-target) | 3.33 | 96.1 | ✓ Strong positive |
| P2 (Jenkins NY1-B04, on-target) | 3.73 | 95.1 | ✓ Strong positive (tighter than cycle-01 full target) |
| P3 (Baker wt1 vs MART-1, mismatch) | 3.55 | 95.4 | Only 0.22 Å degraded vs P1 — consistent with BAKER's own data |
| N1 (scrambled) | 22.0 | 31.4 | ✓ Rejected |
| N2 (random 65aa) | 19.2 | 36.3 | ✓ Rejected |

Dynamic range ~18 Å preserved. Halt thresholds recalibrated to iPAE≤6.0, ipLDDT≥92.0.

### Stage 1 sub-run A results

| Metric | Value |
|---|---|
| Scaffold library | BAKER 152 scaf*.pdb |
| Scaffolds aligned successfully | 152/152 |
| **Bimodal alignment HLA-CA RMSD** | **56 at <1 Å (A\*02:01-native), 96 at 3-6 Å (cross-allele); zero in 1-3 Å middle ground** |
| Designs attempted (partial diffusion, partial_T=15) | 150 |
| Designs completed | 150 (after Trap #33 symlink recovery) |
| **Geometry pass rate** | **72/150 = 48% (3.7× cycle 02's 13%)** |
| Pass population hotspot contacts (median) | **6.5** (cycle 02: ~0) |
| Pass population hotspot contacts (mean) | 7.8 |
| Pass population with ZERO contacts | **0%** (cycle 02: 66%) |
| Fail reasons | 71 insufficient_hotspot_contacts, 15 length_out_of_range |
| **Surprising cross-tab finding (Trap #40)** | **Low-RMSD scaffolds (well-aligned A\*02:01-native) pass at 38%; high-RMSD (cross-allele) pass at 54%. HLA-RMSD is NOT predictive of scaffold transfer quality.** |
| Binder length range (passing) | 66–113 (median 87) |
| Wall time, A100 | 98.5 min |
| Approx GPU cost | ~$5 |
| Partial T schedule | 15 of 50 |
| Motif RMSD during diffusion | ~0.11 Å (excellent preservation) |
| RFdiffusion checkpoint | `Complex_base_ckpt.pt` (cycle 02 sha) |
| Contig per-design format used | `[N-N/0 B1-180/0 C1-9]` |

### Stage 1 sub-run B status

Stubbed (`results/cycle_03/stage1/subrun_b/{subrun_summary.json, designs.jsonl}` empty), deferred. Blocked on: (a) Trap #33 filename fix, (b) sub-run B contig derivation helper.

### Stage 2 truncated status

**Implemented (mock-validated; pending pod real-run).** `splice_binder_subrun_a` composes the 3-chain truncated complex (A=HLA, B=peptide, C=binder); `run_stage2 --target-layout truncated` threads the layout through splice → ProteinMPNN (`--chain_list C`) → FASTA (`HLA:peptide:binder`) → `compute_metrics` (`LAYOUT_CHAINS["truncated"]`). Stage 1 verdict gate accepts sub-run summaries via `_stage1_verdict` (derives PASS from `fraction_geometry_pass >= HALT_THRESHOLD`, imported from `run_stage1`). New Snakemake rule `stage2_subrun_a` (explicit target, not in `rule all`). Real-run command on pod: `snakemake results/cycle_03/stage2/subrun_a/stage2_summary.json --config cycle=03 mock=false -j1`.

### Pull requests in flight

| PR | Branch | Status | Contents |
|---|---|---|---|
| #14 | `claude/cycle-03-prep-g2Eu9` | DRAFT off main | Snakemake DAG, configs, align_scaffolds.py, partial_diffuse.py, contact_filter.py, decomposed iPAE |
| #15 | `feat/cycle-03-baker-truncation` | DRAFT, stacked on #14 | prep_baker_target.py, align_baker_scaffolds.py, LAYOUT_CHAINS truncated, run_controls.py --target=truncated, Trap #31 fix, **Trap #32 contig fix + geometry-gate truncated awareness (latest commit)** |

Neither merged to main. Cycle 03 results live on `feat/cycle-03-baker-truncation`.

### Pod-local hotfixes

| File | Change | Status |
|---|---|---|
| `Snakefile` | `_to_bool()` helper (Trap #34) | ✓ **Committed** (Stage 2 truncated PR) |
| `configs/af2_stage2.yaml` | halt_cut_ipae_max 10.0 → 6.0, halt_cut_iplddt_min 88.0 → 92.0 | ✓ **Committed** (Stage 2 truncated PR) |
| `results/cycle_03/stage1/subrun_a/designs/*.pdb` | 150 symlinks: `design_NNNN.pdb -> design_NNNN_NNNN.pdb` (Trap #33 recovery) | Not committed (gitignored results/); pod state still depends on them for the already-burned run — DO NOT delete. New runs produce single-suffix names natively (Trap #33 fixed), so future symlinks are unnecessary. |

### Recurrence guards needed (test specs for next CC PR)

- ✅ Writer/reader contract test for `partial_diffuse_one` filename output (Trap #33) — **done** (`test_partial_diffuse_one_prefix_has_no_seed`, `test_subrun_a_writer_reader_filename_roundtrip`).
- ✅ Sub-run A real-path Stage 2 mock canary (3-chain design → FASTA → AF2 → metrics) — **done** (`tests/test_stage2_subrun_a.py`).
- Frame-preservation assertion test for `prep_baker_target.py` (preserves coords during truncation).
- Sub-run A real-path Stage 2 mock canary (synthetic 3-chain design → FASTA → AF2 → metrics).

---

*End of state. New thread reading this + the source paper anchors should be able to advise on any technical question without re-deriving prior decisions.*
