# Stage 3 (planned) — ProteinMPNN peptide-context specificity filter

> **Status: PLACEHOLDER — planned for cycle 04, not implemented in the
> cycle-03 prep PR.** This document records the contract so the cycle-03
> architecture has a named landing point for the deferred component.

> **Naming note**: this is distinct from `specs/stage3_crosspan.md` (the in
> silico cross-reactivity panning against off-target peptides). Both touch
> "specificity"; this one is the *upstream design-time* MPNN filter, the
> crosspan spec is the *downstream evaluation*. The "stage3" label here
> follows the kickoff brief's naming and refers to the design-funnel position
> (after the contact filter), not the pipeline Stage 3 of `crosspan`.

## Goal

Bias / filter ProteinMPNN sequence design toward binders that are specific to
the **peptide** surface of the pMHC, not just the conserved MHC framework —
BAKER_LAB_2025's peptide-context specificity step (their `mpnn_spec_filter/`).
This complements cycle-03's:

- amino-acid bias (`configs/proteinmpnn_bias_aa.json`, Trap #30),
- BioPython contact filter (`workflow/scripts/contact_filter.py`),
- decomposed iPAE sub-metrics (`ppi_pae_int_peptide` / `ppi_pae_int_mhc`),

by acting at sequence-design time rather than as a post-hoc geometric or
AF2-based filter.

## Anchored references

- `BAKER_LAB_2025` — `mpnn_spec_filter/`: peptide-context ProteinMPNN scoring
  to enrich peptide-specific binders.
- Cycle 03 — `ppi_pae_int_peptide` / `ppi_pae_int_mhc` decomposition
  (`workflow/scripts/compute_metrics.py`) provides the evaluation-time signal
  this filter should improve at design time.

## Planned inputs / outputs (to be finalized in cycle 04)

- **In**: spliced 4-chain complexes (A/B/C/D), the cycle's MPNN config, the
  BAKER `mpnn_spec_filter/` scoring assets (pod-only).
- **Out**: a per-sequence peptide-specificity score + a filtered
  `sequences.jsonl` feeding the AF2 fan-in.

## Out of scope for cycle 03

No code, configs, fixtures, or rules for this filter ship in the cycle-03 prep
PR. It is recorded here only so the deferred work is tracked and the cycle-03
funnel documents where it will attach (MPNN → spec filter → contact filter →
AF2).
