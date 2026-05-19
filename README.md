# pMHC-I Minibinder Design Pipeline

A reproducible, end-to-end pipeline for de novo minibinder design against
peptide-MHC class I complexes. Mirrors the platform from Johansen et al.
(Science 2025, Hadrup/Jenkins lab) and the Baker lab co-submission
(Liu et al., Science 2025).

**Target:** WT1 RMFPNAPYL / HLA-A*02:01 cancer-testis antigen with broad
clinical interest for AML and solid tumors. Same target as BAKER_LAB_2025,
giving us a direct experimental literature anchor.

**Stack:** RFdiffusion → ProteinMPNN → AF2-multimer (ColabFold) → MD →
in silico cross-panning → active learning. Snakemake-orchestrated;
mock-mode CI; ~40-entry trap book.

**Current state:**
- **Cycle 1 calibration** on 5 published controls established clean
  positive/negative discrimination  ~20 Å iPAE gap and ~60-point ipLDDT
  gap between positives (4.5–4.9 Å, 94–96) and negatives (24.7–25.7 Å,
  31–35). The calibration also reproduced AF2's known peptide-identity
  blind spot (WT1 binder folded well against MART-1), which motivates the
  Stage 3 cross-pan + Hamming proteome scan downstream.
- **Cycle 2 de novo run** on WT1/A*02:01 ran end-to-end: 200 RFdiffusion
  backbones → 26 geometry-pass (13%) → 50 sequences folded by AF2-multimer → 1
  hero design at iPAE 6.41 Å, sitting in the cycle 1 positive-control band.
  2% strict-cut yield matches the pre-optimization range Liu et al. report;
  cycle 03 applies peptide-centric arcing, alanine bias, and increased sequence
  sampling per backbone for the production run.

📖 [Scientific narrative](docs/narrative.md)
📊 [Cycle 02 results](notebooks/cycle_02_report.ipynb)
🪤 [Known traps](docs/known_traps.md)
⚙️ [Reproducibility](docs/narrative.md#engineering-rigor-and-reproducibility)
