# Stage 2 — ProteinMPNN sequence design + AF2 fan-in (real designs)

> **Controls pipeline (P1–N2 literature binders, cycle 1 contract) is a
> separate spec: see [`specs/stage2_controls.md`](stage2_controls.md).** The
> controls panel runs on every cycle and provides cross-pipeline calibration;
> this file covers Stage 1 outputs only.

## Goal

Take Stage 1's 200 RFdiffusion backbones, design 4 sequences per backbone
with ProteinMPNN (800 sequences total), splice each as chain D onto the
cleaned pMHC (A=HC, B=β2m, C=peptide), rank by MPNN NLL, fold the top 100
with AF2-multimer at `num_recycles=6`, score iPAE/ipLDDT/BSA, apply the
halt gate.

## Anchored references

- `HADRUP_JENKINS_2025` — Fig 1A pipeline: 5500 backbones × 4 seqs → AF2 →
  intermediate cut `iPAE<12 ipLDDT>88` → 44 designs at the tight cut
  `iPAE<7 ipLDDT>92`. Stage 2 mirrors the first half of this funnel.
- `BAKER_LAB_2025` — Methods: 4–32 MPNN seqs/backbone; ProteinMPNN with
  peptide-context for specificity (deferred to our Stage 3).
- Cycle 1 — Stage 2 controls validated (positives iPAE 4.5–4.9 Å, negatives
  24.7–25.7 Å). Same AF2 scoring pipeline inherited.

## Inputs

1. `results/cycle_NN/stage1/rfdiffusion/stage1_summary.json` — verdict must
   be `PASS` (orchestrator exits 1 otherwise).
2. `results/cycle_NN/stage1/rfdiffusion/designs.jsonl` — 200 per-design records.
3. `results/cycle_NN/stage1/rfdiffusion/designs/design_*.pdb` — backbones,
   binder = chain "A" internally.
4. `results/cycle_NN/stage0/target.yaml` — TargetManifest. Splicing uses
   `target.primary.cleaned_pdb` as the chain-A/B/C source.
5. `configs/proteinmpnn_default.yaml`:
   ```yaml
   model: v_48_020
   sampling_temperature: 0.1
   n_seqs_per_backbone: 4
   fixed_chains: [A, B, C]
   design_chains: [D]
   ```
6. `configs/af2_stage2.yaml`:
```yaml
   # Cycle 02 PoC overrides — see "Cycle 02 PoC scaling" section below.
   num_recycles: 3              # cycle 02 PoC; cycle 03+ promotes top survivors to 6
   fan_in_top_n: 50             # cycle 02 PoC; raise to 100 once Stage 2 is validated
   halt_cut_ipae_max: 12.0
   halt_cut_iplddt_min: 88.0
   halt_threshold_fraction: 0.10
```
7. `configs/seeds.yaml` — `formulas.proteinmpnn` and `formulas.af2_fanin`
   keys are canonical; `reserved.cycle_02_proteinmpnn` and
   `reserved.cycle_02_af2_fanin` are assertion bounds.
8. `/workspace/ProteinMPNN/` (bootstrap.sh clones).
9. `/workspace/ProteinMPNN/vanilla_model_weights/v_48_020.pt` (bundled in
   the clone; bootstrap.sh sha256-verifies in place).

## Outputs

```
results/cycle_NN/stage2/
├── proteinmpnn_designs/
│   ├── sequences.jsonl            # 800 records: backbone_id, seq_id (0..3), seq, mpnn_nll, mpnn_seed, backbone_pdb, spliced_pdb
│   └── spliced/
│       └── design_NNNNN.pdb       # 200 4-chain spliced complexes (A=HC, B=β2m, C=peptide, D=binder)
├── af2_designs/
│   ├── predictions/
│   │   └── design_NNNNN_seqK/     # ColabFold per-prediction subdir (ranked_0.pdb + scores.json)
│   └── metrics.jsonl              # 100 records: design_id, seq_id, iPAE, ipLDDT, BSA, plddt_chain, pae_chain
├── stage2_summary.json            # halt verdict, counts, top-10 ranking
└── run.log
```

`stage2_summary.json` schema:

```json
{
  "cycle": 2,
  "target_id": "wt1_a0201",
  "n_backbones_in": 200,
  "n_sequences_designed": 800,
  "n_af2_folded": 50,
  "n_pass_intermediate": 8,
  "fraction_pass_intermediate": 0.16,
  "fraction_pass_tight": 0.04,
  "halt_rule": {
    "name": "fraction_pass_intermediate",
    "threshold": 0.10,
    "observed": 0.17,
    "verdict": "PASS",
    "note": "Threshold is calibration-only for cycle 2; tighten in cycle 3"
  },
  "top_10_by_ipae": [
    {"design_id": "design_00042", "seq_id": 1, "iPAE": 5.3, "ipLDDT": 91.2, "BSA": 1320.5}
  ],
  "wall_minutes_total": 540,
  "mock": false
}
```

Both `fraction_pass_intermediate` (Jenkins intermediate cut) and
`fraction_pass_tight` (Jenkins final cut `iPAE<7 ipLDDT>92`) are reported.
Only the intermediate fraction drives the halt; tight is for cycle 3
calibration.

## Quality gate (halt rule)

**Rule:** `fraction_pass_intermediate >= 0.10`, calibration-only for cycle 2.

Per-design pass condition: `iPAE < 12.0 AND ipLDDT > 88.0` (Jenkins intermediate).

Threshold rationale: Jenkins's end-to-end yield is <0.2 % at tight cuts
(44/22000 sequences); Householder ~4 % at the scaffold level. At our looser
intermediate cuts on a 100-design AF2 fan-in, **0.10 is a pipeline-health
calibration** — captures "the pipeline produces plausibly-foldable
designs", not "the pipeline produces winners". If `observed < 0.10`,
something upstream is broken (Stage 1 hotspots wrong, MPNN config wrong,
chain splicing wrong) — halt loud. Tighten in cycle 3.

The halt gate uses `>=` (PASS at boundary). The mock fixture deliberately
constructs 4/40 pass exactly to exercise this boundary; the corresponding
test in `tests/test_stage2_designs_halt_gate.py` is a refactor canary.

## Cycle 02 PoC scaling

For the first end-to-end real Stage 2 run (cycle 02), we reduce two
parameters from the asymptotic spec defaults to keep wall time tractable:

| Parameter         | Asymptotic | Cycle 02 PoC | Rationale |
|-------------------|------------|--------------|-----------|
| `num_recycles`    | 6          | **3**        | BindCraft / dl_binder_design default for binder triage. Halves AF2 wall time. Cycle 03+ promotes top survivors to 6 recycles for sharper final metrics. |
| `fan_in_top_n`    | 100        | **50**       | At expected 10 % intermediate pass rate this yields ~5 passing designs — enough signal to validate the pipeline produces real binders, accepting wider statistical uncertainty on the halt verdict (false-PASS rate ~24 % at true 5 % rate). Cycle 03 raises to 100 for a confident gate verdict. |

Total Stage 2 wall time at cycle 02 settings: **~3.5 h** end-to-end
(~30 min ProteinMPNN + ~2.5 h AF2 + overhead) vs ~11 h at asymptotic
defaults.

If cycle 02 PASSES halt at fan-in 50 and you want a firmer verdict before
committing to cycle 03 architecture changes, re-run AF2 on the remaining
50 predictions (designs ranked 51–100 by MPNN global_score) — this is a
cheap follow-on, not a redo.
## Mock mode

`uv run python scripts/run_stage2.py --mock --cycle 99` reads Stage 1's
mock fixtures (10 backbones), generates 4 deterministic mock MPNN sequences
per backbone, splices, "folds" via mock AF2 (copies pre-generated mock
prediction outputs from `tests/fixtures/stage2/designs/`), and emits
`stage2_summary.json`. Fixture calibration: exactly **4 of 40 mock AF2
predictions pass** the intermediate cut → `observed == 0.10` (boundary,
PASS by `>=`).

The fan-in top-N is config-driven: real mode uses
`af2_stage2.yaml.fan_in_top_n=100`, mock mode caps at `min(40, n_seqs)`
so the boundary fixture works without scaling up.

`pytest` additionally exercises a synthetic fail-path: 3/40 passing →
`observed == 0.075` → verdict `FAIL`.

## Chain ordering contract (inherited from Stage 1)

A=HC, B=β2m, C=peptide, D=binder. Stage 1 emits binder as chain "A"
internally. `workflow/scripts/splice_binder.py` renames to "D" and writes
a 4-chain PDB; `scripts/run_stage2.py` re-asserts the chain ordering at the
top of the AF2 input FASTA writer. **Never pass Stage 1 PDBs directly to
ColabFold** — they'd be silently treated as 1-chain HC inputs. See
`docs/known_traps.md` trap #17.

## Seed reservations (configs/seeds.yaml)

| Stage         | Formula                                            | Cycle 2 range |
|---------------|----------------------------------------------------|---------------|
| ProteinMPNN   | `cycle * 1000 + 200 + design_index * 4 + seq_index`| 2200–2999     |
| AF2 fan-in    | `cycle * 1000 + 1000 + fan_in_rank`                | 3000–3099     |

`workflow/scripts/aggregate_mpnn_outputs.py` exposes `_compute_mpnn_seed`
and `_compute_af2_seed`, both of which assert range membership against
`configs/seeds.yaml::reserved`. Tests in
`tests/test_aggregate_mpnn_outputs.py` cover range + uniqueness across all
800 MPNN records and all 100 AF2 records.

## Known traps

1. **Chain renaming** (trap #17): Stage 1 emits binder as chain "A".
   Stage 2 must rename to "D" before AF2. The splice helper handles this;
   never skip it.
2. **ProteinMPNN's `--chain_id_jsonl` and `--fixed_positions_jsonl`** are
   JSON-encoded but format-finicky. Mirror the cycle 1 controls' usage in
   `workflow/scripts/run_proteinmpnn.py`; consult `dauparas/ProteinMPNN`
   examples if extending.
3. **Cleaned PDB residue ranges must match Stage 1's contigmap exactly.**
   If Stage 1 was run against `3hpj_clean.pdb` but Stage 2 reads a
   regenerated cleaned PDB with different residue numbering, splicing
   fails. The splice helper always reads `target.primary.cleaned_pdb`
   from the cycle's `stage0/target.yaml` — same source-of-truth pattern
   as Stage 1.
   4. **`--use-gpu-relax` skipped during AF2 triage.** Baker lab's
   `pmhc_fold.py` uses `do_relax=False` during high-throughput screening.
   Drop the flag from the `colabfold_batch` invocation; saves ~30–50 % wall
   time. Apply AMBER relax only to top-N final structures in cycle 03+ if
   BSA quality becomes a concern.
5. **ProteinMPNN `.fa` output starts with the original (input) sequence.**
   Aggregator must skip the first record per `.fa` file or you'll get
   `n_seqs_per_backbone + 1` records instead of `n_seqs_per_backbone`.
6. **ProteinMPNN outputs the full complex sequence** (all four chains) per
   record, in PDB chain order. The aggregator must slice the binder
   portion using the known binder length from the Stage 1 manifest; chains
   A/B/C in the output match the native sequence and are discarded.
7. **ColabFold multimer FASTA is colon-separated within a single record**,
   not multi-record one-chain-per-entry. The multi-record format is
   DeepMind's original AF2-Multimer convention; `colabfold_batch` requires
   colons: `>id\n<chainA>:<chainB>:<chainC>:<chainD>\n`. No trailing colon.

## Done criteria

- [ ] `uv run python scripts/run_stage2.py --mock --cycle 99` exits 0 in <30 s, verdict `PASS`.
- [ ] `pytest tests/test_stage2_designs_halt_gate.py` green (pass + fail paths).
- [ ] No absolute paths in `tests/fixtures/stage2/designs/*` (recurrence guard).
- [ ] CI green on `uv sync --extra dev`.
- [ ] Fresh-clone preflight green before push.
- [ ] Real pod run on cycle 02 deferred to post-merge; this PR ships mock-CI-green only.

---

## Cycle 03 changes

Cycle 03 keeps the MPNN → splice → AF2 fan-in shape but tightens it and adds a
placement gate. All changes are real-mode (the mock 4/40 halt-gate calibration
is preserved unchanged).

### ProteinMPNN amino-acid bias (`configs/proteinmpnn_cycle03.yaml`)

RFdiffusion's crude-sequence prior is Ala-heavy (Trap #30). Cycle 03 selects
`proteinmpnn_cycle03.yaml` (via `MPNN_DESIGNS_CONFIG` in the Snakefile, cycle
≥ 3) which adds `bias_AA_jsonl: configs/proteinmpnn_bias_aa.json`
(`A: -2.0`, `E/L/R: +1.0`, BAKER redesigned-binder composition reference).
`scripts/run_stage2.py::_run_real_proteinmpnn` passes `--bias_AA_jsonl` when
the config carries it.

### Contact filter (`workflow/scripts/contact_filter.py`)

A Rosetta-free BioPython gate inserted **between MPNN and AF2** in the real
funnel: keep a design only if its binder (chain D) makes ≥ `min_contacts`
C-beta atoms within `distance` Å (default 5.0 Å) of any peptide (chain C)
atom. This discards mis-placed binders cheaply, before the expensive AF2
fan-in — the direct lever against the cycle-02 placement deficit. It is
unit-tested and spec'd here; it is intentionally **not** inserted into the
mock Snakemake DAG so the existing Stage 2 4/40 halt calibration is untouched.

### Decomposed iPAE (`workflow/scripts/compute_metrics.py`)

BAKER's `ppi_pae_int` is our combined iPAE (D ↔ A+B+C). `metrics_for_design`
now also emits `ppi_pae_int_peptide` (D ↔ C) and `ppi_pae_int_mhc` (D ↔ A+B),
so we can tell binders that grip the specificity-bearing peptide apart from
binders that only hold the conserved MHC framework. Surfaced in
`af2_designs/metrics.jsonl` (real mode).

### Tighter funnel (`configs/af2_stage2.yaml`)

Promoted from the cycle-02 PoC to the asymptotic spec: `num_recycles = 6`,
`fan_in_top_n = 100`, and the intermediate halt cut tightened to
`halt_cut_ipae_max = 10.0` (was 12.0). The halt rule and `>=` boundary
semantics are otherwise unchanged.

### Deferred to cycle 04

The ProteinMPNN peptide-context **specificity filter** (`mpnn_spec_filter/`)
is planned but not implemented here — see `specs/stage3_spec_filter.md`.
