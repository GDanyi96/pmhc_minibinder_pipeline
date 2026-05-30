# pMHC-I Minibinder Design Pipeline

A reproducible de novo minibinder design against
peptide-MHC class I complexes, with a peptide-resolved interface metric that
measures whether a design reads the disease-defining peptide or only the
conserved MHC framework. 

**Target:** WT1 RMFPNAPYL / HLA-A\*02:01, a cancer-testis antigen with broad
clinical interest for AML and solid tumors. Same target as the Baker lab paper,
giving a direct experimental literature anchor.

![figures/design_3010_peptide_interface.png](docs/figures/design_3010_peptide_interface.png).

**Stack:** RFdiffusion → ProteinMPNN → AF2-multimer (ColabFold) → interface-metric
decomposition → (planned) MD, in-silico cross-panning, active learning.
Snakemake-orchestrated; mock-mode CI; all thresholds externalized to config.

## Current state

- **Cycle 1 — calibration.** Five published controls established clean
  positive/negative discrimination on the full-target assay (positives iPAE
  4.5–4.9 Å / ipLDDT 94–96; negatives 24.7–25.7 Å / 31–35). The deliberate
  specificity-failure control (a WT1 binder folded against MART-1) reproduced
  AF2's known peptide-identity blind spot on aggregate metrics — but degrades on
  the peptide-restricted channel, motivating the metric below.
- **Cycle 2 — de novo run.** End to end on WT1/A\*02:01: 200 RFdiffusion
  backbones → 26 geometry-pass (13%) → 50 AF2 folds → 1 design crossing the
  strict cuts. The peptide-contact decomposition reclassified that "hero" as a
  **peptide-blind MHC-framework binder** (closest residue 28–40 Å from the
  peptide). Learning used for next cycles.
- **Cycle 3 — BAKER-scaffold partial diffusion.** 150 scaffolds → 72
  geometry-pass (48%) → 288 sequences → 288 AF2 folds. Only 32% engage the
  peptide; **one** design (design_3010_seq00) reads both WT1 specificity residues
  (N5, Y8) in the control-grade band, via a single binder residue (R55).
  Framework bias is structural (a charge hypothesis was raised and falsified),
  so the cycle-04 fix is at the generation stage.

  



All results are in silico predictions; no experimental validation has been
performed.

📖 [Scientific narrative](docs/narrative.md)
🧪 [Cycle 02 results (de novo)](docs/cycle_02.md)
🧪 [Cycle 03 results (BAKER scaffolds)](docs/cycle_03.md)
🎯 [3010 peptide interface (figure)](docs/figures/design_3010_peptide_interface.png)
🧠 [Methodological lessons](docs/methodological_lessons.md)
📊 [Analysis tables](results/cycle_03/analysis/)
🔧 [Engineering notes (reproducibility)](docs/known_traps.md)
