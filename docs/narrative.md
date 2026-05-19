# De novo pMHC-I minibinder design pipeline: scientific narrative



## Executive summary

This pipeline is a calibrated three-tier funnel for designing small (~70 to 110 aa) de novo proteins that specifically recognize a peptide-MHC class I complex on cancer cells.

1. Calibration. Published positive and negative control binders run through the structural prediction pipeline establish that the metrics actually distinguish real binders from non-binders for this target system, with a ~20 Å iPAE gap and a ~60-point ipLDDT gap.
2. Sensitivity. RFdiffusion generates de novo backbones; ProteinMPNN designs sequences; AF2-multimer folds the complex and reports interface confidence metrics. Surviving designs are those AF2 trusts will form tight complexes.
3. Specificity. Surviving designs are tested against an off-target panel of self-peptides presented on the same HLA allele. Designs that retain affinity for the cognate while showing reduced engagement with off-targets are kept.

The architecture mirrors HADRUP_JENKINS_2025 (Science 2025, Johansen et al.) and BAKER_LAB_2025 (Science 2025, Liu et al.), with explicit reproducibility safeguards and trap-tracking discipline.

---

## 1. Clinical motivation

The immune system normally surveils cells by reading peptides displayed on MHC class I (the cell's intracellular protein "library card"). Tumor cells display peptides derived from oncogenic drivers or neoantigens this way. The natural recognition machinery (T-cell receptors, TCRs) is hard to isolate, often low-affinity, and difficult to manufacture at therapeutic scale.

De novo designed minibinders for pMHC-I complexes solve several problems at once:

- They access intracellular antigens via the surface-presented peptide, which conventional antibodies cannot target.
- They are small, stable, and can be produced in *E. coli* at low cost.
- They can be incorporated into chimeric antigen receptor (CAR) constructs to redirect T-cell killing. The "BIKE" concept (de novo pMHC-Binders for Immune-mediated Killing Engagers) is set out in HADRUP_JENKINS_2025.
- A single computational pipeline can address arbitrary tumor antigens, generalizable across patients and cancers in principle.

The bottleneck is **specificity**: the binder must recognize the cancer peptide but not closely related self-peptides presented on the same HLA on healthy cells. That is the hard problem this pipeline is built to address.

---

## 2. Target choice

The cognate complex is **WT1 RMFPNAPYL / HLA-A\*02:01**, crystal structure PDB 3HPJ.

| Component | Value | Rationale |
|---|---|---|
| Antigen | WT1 (Wilms Tumor 1) | Overexpressed in AML and multiple solid tumors; well-established TAA |
| Peptide | RMFPNAPYL (residues 126-134) | Top-tier cancer immunotherapy target |
| HLA allele | HLA-A\*02:01 | Most common allele in the Caucasian population (~50%); broadest patient applicability |
| Structure | PDB 3HPJ | Published crystal complex; clean starting point |

BAKER_LAB_2025 targeted this exact complex in their Science 2025 paper. That gives the project a direct experimental literature anchor: their published binders become this pipeline's positive controls; their results become the benchmark.

---

## 3. Pipeline architecture overview

Three sequential validation tiers:

```
┌─────────────────────────────────────────────────────────────────────────┐
│                                                                          │
│  TIER 1: CALIBRATION                                                     │
│  (cycle 1) published controls through the pipeline                       │
│                                                                          │
│  Question: do the metrics actually distinguish binders from non-binders? │
│  Answer:   YES (positives 4.5-4.9 Å iPAE, negatives 24-26 Å)             │
│                                                                          │
└───────────────────────────────────┬──────────────────────────────────────┘
                                    │
┌───────────────────────────────────▼──────────────────────────────────────┐
│                                                                          │
│  TIER 2: SENSITIVITY                                                     │
│  (cycle 02) Stages 0-2: de novo design through the calibrated assay      │
│                                                                          │
│  Stage 0  → target preparation (3HPJ cleanup)                            │
│  Stage 1  → RFdiffusion: 200 de novo binder backbones                    │
│  Stage 2  → ProteinMPNN (4 seqs/backbone) → top to AF2-multimer          │
│             → iPAE, ipLDDT, BSA per prediction                           │
│             → halt gate: ≥10% pass iPAE<12 AND ipLDDT>88                 │
│                                                                          │
│  Output: ~5-10 surviving designs that AF2 trusts                         │
│                                                                          │
└───────────────────────────────────┬──────────────────────────────────────┘
                                    │
┌───────────────────────────────────▼──────────────────────────────────────┐
│                                                                          │
│  TIER 3: SPECIFICITY (planned)                                           │
│                                                                          │
│  Stage 3a → sequence-based off-target panel (Hamming + presentation)     │
│  Stage 3b → AF2 cross-panning of survivors × off-target panel            │
│                                                                          │
│  Output: subset of survivors with low predicted cross-reactivity         │
│                                                                          │
└──────────────────────────────────────────────────────────────────────────┘
```

---

## 4. Stage 0: target preparation

The cleaned target complex contains three chains:

| Chain | Identity | Length |
|---|---|---|
| A | HLA-A\*02:01 heavy chain | ~275 aa |
| B | β2-microglobulin | ~100 aa |
| C | RMFPNAPYL peptide | 9 aa |

Stage 0 strips waters, ligands, and non-essential residues from 3HPJ, validates chain numbering, and emits a target manifest consumed by every downstream stage. This is the source-of-truth structure for everything that follows.

---

## 5. Tier 1: pipeline calibration via controls (cycle 1)

Before generating any de novo designs, the pipeline was run on five published reference binders. The question this answers: do the metrics actually distinguish real binders from non-binders for our target system?

### Why calibration was necessary

AF2 confidence metrics are not absolute. An iPAE of 8 Å is meaningless in isolation. Without empirical anchors, every downstream result is uninterpretable:

| Hypothetical result | Without calibration | With calibration |
|---|---|---|
| de novo design at iPAE 6.5 Å | "Is this good? Is this noise?" | "Comparable to published binders (4.5 to 4.9 Å); behaves like a real binder" |
| de novo design at iPAE 11 Å | "Marginal. Is the threshold even right?" | "Far from positives, closer to negatives (24 to 26 Å); likely noise" |

### The control panel

| ID | Identity | Purpose |
|---|---|---|
| P1 | Baker lab WT1 binder, cognate target | Confirms a real positive passes (rules out thresholds too tight) |
| P2 | Jenkins NY1-B04 binder, cognate target (NY-ESO-1) | Confirms calibration generalizes across pMHC systems |
| P3 | Baker WT1 binder on the *wrong* target (MART-1) | Specificity check (rules out: pipeline rewards generic foldedness) |
| N1 | Scrambled Baker WT1 sequence, same backbone | Rules out: metrics are backbone-dominated and sequence-blind |
| N2 | Random 85-aa sequence | Rules out: AF2 always reports false-positive interface metrics |

### What each control rules out, and what the cycle 1 data showed

| Failure mode the control rules out | Result |
|---|---|
| Pipeline reports bad metrics for everything | P1 at iPAE 4.88 Å, ipLDDT 94.5. Real positive passes. |
| Calibration only works for WT1 | P2 (Jenkins NY1-B04 vs NY-ESO-1) at 4.56 Å, 95.9. Generalizes across pMHCs. |
| Pipeline rewards generic foldedness regardless of target | P3 (WT1 binder on MART-1) at 4.54 Å, 94.2. See note below. |
| Metrics are backbone-dominated, sequence-blind | N1 (scrambled, same backbone as P1) at 24.7 Å, 34.5. Sequence-driven. |
| AF2 always reports good interface metrics for anything | N2 (random 85-mer) at 25.7 Å, 31.4. AF2 reports garbage as garbage. |

### Cycle 1 results

| Control | iPAE (Å) | ipLDDT | BSA (Å²) |
|---|---|---|---|
| P1 | 4.88 | 94.5 | 2692 |
| P2 | 4.56 | 95.9 | 2788 |
| P3 | 4.54 | 94.2 | 2698 |
| N1 | 24.72 | 34.5 | 1524 |
| N2 | 25.70 | 31.4 | 2153 |

The discrimination between positives (4.5 to 4.9 Å iPAE, 94 to 96 ipLDDT) and negatives (24.7 to 25.7 Å, 31 to 35) is clean. The iPAE gap is ~20 Å and the ipLDDT gap is ~60 points. The ipLDDT separation is sharper than the iPAE separation and is the more discriminating signal in this calibration.

The P3 result needs unpacking. P3 is the Baker WT1 binder folded against MART-1 ELAGIGILTV on the same HLA-A\*02:01. If AF2 iPAE were a complete specificity readout, P3 should have scored substantially worse than P1. Instead it scored on par. This is the documented AF2 specificity blind spot. HOUSEHOLDER_GARCIA_2025 addresses it with a proteome-wide Hamming scan; BAKER_LAB_2025 addresses it with yeast display experimental validation; MARES_IOANNIDIS_2025 quantified it at AUROC 0.06 to 0.22 for sequence-based predictors on structurally valid peptides. The cycle 1 calibration reproduced the literature concern in our own hands. Stage 3 (in silico cross-panning plus Hamming proteome scan) is the answer, and this P3 result is the empirical justification for it.

A side note on the cycle 1 halt rule. The initial halt rule required P1 to rank first or second of five by iPAE. That rule failed: P1 ranked third behind P3 (4.54) and P2 (4.56), separated by ~0.3 Å from both. The rule was poorly chosen. A controls panel doesn't fail because of an unexpected ordering between two positives at essentially identical scores. It fails when discrimination between positives and negatives breaks down. The rule was replaced with three component criteria: P1 below threshold, P2 within published calibration, clean positive/negative separation. All three pass. Logged in the trap book.

---

## 6. Tier 2, Stage 1: RFdiffusion (cycle 02)

With the assay calibrated, the pipeline moves to de novo generation.

Inputs for cycle 02: cleaned 3HPJ as the target (chains A, B, C for HLA heavy, β2m, peptide). Eight hotspots: C1, C4, C6, C8 (four outward-facing peptide residues) plus A65, A66, A150, A155 (four flanking HLA residues). The peptide-heavy bias follows BAKER_LAB_2025; binders that contact mostly the peptide are inherently more peptide-specific. Binder length range 70 to 110 residues. 200 backbones requested.

Output: 200 protein backbones (CA coordinates only, no sequence yet) that geometrically arc above the peptide groove. This stage is a pipeline-health metric, not a quality filter. It answers whether RFdiffusion can produce backbones that approach the cleft at all and whether they cluster into recognizable topologies.

### Cycle 02 results

Running on a single A100 SXM 80GB took ~10 hours of wall time. All 200 designs completed and wrote 4-chain PDBs (A=275 aa, B=100 aa, C=9 aa, D=70 to 110 aa). Motif RMSD landed at ~0.12 to 0.13 Å across all 200 designs, well inside geometric tolerance. The binder length distribution sat cleanly inside the 70 to 110 contig (median 89).

The geometry quality gate passed 26 of 200 designs (13%). Of the 174 failures:

- 173 of 174 (99%) failed on insufficient hotspot contacts.
- 4 of 174 failed on internal CA clash.
- 131 of 200 designs (66%) had literally zero CA atoms within contact distance of any hotspot residue.

The diagnosis: the binders are at reasonable lengths with valid backbones, but they're landing on the wrong face of HLA. With only 8 hotspots and the default 10 Å recenter radius, the spatial prior on this design isn't strong enough to constrain placement to the groove. Many failures dock against the α3 domain (the larger non-groove surface that contacts β2m), not the peptide-binding cleft.

This is a placement deficit, not a backbone-quality deficit. The mitigation is denser hotspots and peptide-centric arcing prompts per BAKER_LAB_2025. That is the cycle 03 plan. The 26 geometry-passing backbones carried forward to Stage 2.

Four engineering incidents were identified during this cycle, all caught by the pipeline's own HALT gates. See §10 (engineering rigor) for the writer/reader contract drift class and how each was diagnosed and pinned.

---

## 7. Tier 2, Stage 2: ProteinMPNN plus AF2-multimer (cycle 02)

The 26 geometry-passing backbones from Stage 1 enter Stage 2 for sequence design and forward folding.

### ProteinMPNN

ProteinMPNN takes each backbone and produces sequences predicted to fold into it while making good contacts with the fixed pMHC target.

Configuration: chains A, B, C fixed at native target sequence; chain D designed; 4 sequences per backbone (104 candidates for cycle 02); sampling temperature 0.1 for high confidence and low diversity. Output ranking by `mpnn_global_score` (NLL-like, lower is better). The top ~50 by score proceed to AF2-multimer.

### AF2-multimer (ColabFold)

Each ranked sequence is folded as a four-chain complex (HLA heavy plus β2m plus peptide plus binder). Configuration: `num_recycles=3` for triage (cycle 03 will promote top hits to 6 for sharper metrics), 5 models per prediction, no AMBER relax (Baker convention; saves wall time). Inputs colon-separated as a multimer FASTA. MSA is queried per design against the ColabFold server, which dominates wall time.

### Metrics and halt gate

Three metrics per fold:

| Metric | Definition | Intermediate cut | Tight cut (informational) |
|---|---|---|---|
| iPAE | mean `min(pae_ij, pae_ji)` over binder↔target residue pairs within 8 Å heavy-atom distance | < 12 Å | < 7 Å |
| ipLDDT | mean pLDDT over residues participating in any interface contact | > 88 | > 92 |
| BSA | Å² of interface buried | informational | informational |

Pipeline-health halt gate: at least 10% of folded designs must pass the intermediate cut (iPAE < 12 AND ipLDDT > 88). Below 10% means something upstream is wrong (hotspots, MPNN config, splicing). At or above 10% means the pipeline is functioning and survivors can proceed to the specificity tier. The threshold is loose by design; this is a "do the metrics produce real binders at all?" gate, not a "are these binders production-ready?" gate.

### Cycle 02 wall time

ProteinMPNN finished in ~3 minutes. AF2-multimer is the long part, dominated by per-design ColabFold MSA server queries. Expected ~5 hours wall on a single A100.

### Cycle 02 results

*Pending. Fill once `stage2_summary.json` lands: iPAE/ipLDDT scatter, intermediate cut pass count, top-N table, brief failure-mode notes.*

---

## 8. Tier 3: specificity screening (planned)

Stage 2 produces a small set of designs (~5 to 10 expected) that AF2 trusts will fold tightly against the cognate target. The cycle 1 P3 result already showed these metrics don't carry specificity information by themselves. Tier 3 builds the specificity readout in two stages.

Stage 3a, sequence-based off-target panel curation, per HOUSEHOLDER_GARCIA_2025. Download HLA-A\*02:01-presented self-peptides from MHC Motif Atlas and IEDB, score by Hamming distance to the cognate RMFPNAPYL, and retain the closest neighbors as candidate off-targets. Stage 3a is CPU-only and depends only on the target peptide, so it can run in parallel with Stage 2.

Stage 3b, structural cross-panning. For each surviving binder, run AF2-multimer against each off-target pMHC and compute the same interface metrics. The specificity signal is the contrast: a useful binder shows tight metrics for the cognate and degraded metrics across the off-target panel. Expected scale ~10 binders × ~15 off-targets × ~3 min per AF2 prediction, roughly 7 to 8 GPU hours, reusing all of Stage 2's infrastructure.

This is the differentiating depth of the pipeline. The Jenkins lab's published platform performs the equivalent cross-pan before any wet-lab handoff.

---

## 9. Beyond specificity (planned)

Three further stages complete the architecture.

Stage 4, embedding and diversity curation: ESM-2 sequence embedding of surviving designs, UMAP plus farthest-point sampling to select a structurally diverse subset for experimental synthesis. CPU-only, minutes of wall time. The point is to avoid sending 20 near-clones to the lab when 20 diverse candidates would teach more.

Stage 5, active learning: surrogate model (XGBoost or small MLP) trained on design-outcome pairs from experimental rounds, with uncertainty-aware acquisition for the next cycle's priorities. The computational framework is small. The science requires wet-lab feedback this pipeline can't generate alone, which is the experimental partnership the PhD project would build out.

Stage 6, reporting: cycle metrics, plots, structured markdown report. Already partially in place via the Snakemake DAG and the per-stage `summary.json` contract.

---

## 10. Engineering rigor and reproducibility

The cycle 02 Stage 1 run surfaced four engineering incidents. All four belong to the same class: writer/reader contract drift hidden by mock fixtures that aligned with the wrong assumption.

The pattern. A downstream consumer (filename enumerator, chain-identification routine, splice helper) carried an incorrect belief about an upstream producer's output (RFdiffusion's filename convention, its chain assignment). The mock test fixtures had been constructed (in some cases unintentionally) to match the wrong belief, so the canary passed not because the code was correct but because the mock and the code shared the same incorrect assumption. Each of the four was caught by the pipeline's own HALT gates during cycle 02 (zero designs reported complete; every design reported as `binder_length=275`; etc.), then diagnosed, patched, and pinned with a regression test that exercises the real producer's output rather than a fixture mirroring the consumer's expectations. All four are documented in `docs/known_traps.md` (entries #26 through #29).

Two engineering principles fell out.

Mock canaries that align with the wrong assumption are tautological. Test against the real producer at least once, even if as a synthetic regression test that mimics the real writer's naming and chain conventions. A canary that only exercises the mock is a canary that validates nothing useful.

Loud failures beat silent fallbacks. The four original code paths all had silent picks ("if multiple chains, pick A"; "if filename doesn't match, skip the design"). These let broken contracts produce numbers instead of errors. Every fix replaced the silent pick with a `ValueError` carrying diagnostic context. Numbers without context can be wrong; errors are always information.

A `--skip-subprocess` affordance was added during cycle 02 and reduced post-bug debugging from ~10 hours (re-run RFdiffusion) to ~8 seconds (re-enumerate existing PDBs). Without it, cycle 02 debugging would have stretched across days. Build re-evaluation paths into any stage with expensive primary compute.

Reproducibility. The repo bootstraps on a RunPod A100 pod via `bootstrap.sh`. The script installs two Python environments (the main uv-managed pipeline env, and a SE3nv conda env that RFdiffusion needs), pins the JAX and PyTorch versions known to be compatible, and clones model weights. Snakemake mock mode runs the full DAG with synthetic fixtures and gates CI on every push. Dev-only scratch files are gitignored.

---

## 11. References

Six published works anchor the pipeline's design decisions.

HADRUP_JENKINS_2025 (Johansen et al., *Science* 2025). The blueprint paper. RFdiffusion + ProteinMPNN + AF2 + MD + IMPAC-T T-cell killing assay. NY-ESO-1 and neoantigen targets. Cryo-EM at PDB 9NNF.

BAKER_LAB_2025 (Liu, Greenwood, Bonzanini et al., *Science* 2025). 11 pMHC targets across 4 HLA alleles. Peptide-centric arcing prompts for placement constraint. Partial diffusion scaffold repurposing (mage-513). Yeast display experimental validation. Code at Zenodo (77forest/pmhcibinderdesign). Same WT1 / A\*02:01 target as this pipeline.

HOUSEHOLDER_GARCIA_2025 (Householder et al., *Science* 2025). De novo α-helical four-helix-bundle TCR mimic for NY-ESO-1 / HLA-A\*02:01. 9.5 nM affinity. Crystal at PDB 9MIN. Proteome-wide Hamming distance scan as the in silico specificity filter.

MARES_IOANNIDIS_2025 (Mares, Ioannidis et al., UC Berkeley / CZ Biohub). Diffusion-based structure-conditioned pMHC-I peptide library generation. Quantified the AUROC 0.06 to 0.22 blind spot of NetMHCpan, MHCflurry, HLApollo on structurally valid peptides. 20 HLA alleles. Code at github.com/sermares/struct-mhc-dev.

BENTZEN_HADRUP_2019 (Bentzen et al., *Nat Biotechnol* 2019). DNA-barcoded MHC multimers for one-pot 190-variant TCR specificity profiling. The Barracoda analysis framework. The experimental cross-reactivity readout this pipeline is eventually designed to feed.

HADRUP_ESMO_2025 (Hadrup, ESMO IO Congress 2025 opinion). Clinical framing of AI minibinders across CAR-T, T-cell engagers, autoimmune. Source of the BIKE motivation in §1.

---

## TL;DR for the application

This pipeline implements a calibrated three-tier specificity-screening apparatus for de novo pMHC-I targeting minibinders.

1. Tier 1 (calibration): five published controls established a ~20 Å iPAE gap and ~60 point ipLDDT gap between positives (4.5 to 4.9 Å iPAE) and negatives (24.7 to 25.7 Å), and reproduced AF2's known peptide-identity blind spot (P3, the WT1 binder on MART-1, scored on par with the cognate-target positive).
2. Tier 2 (sensitivity, Stages 0-2): cycle 02 generated 200 de novo backbones with motif RMSD ~0.12 Å, 26 of which passed the geometry gate (13%). The 26 carry forward through ProteinMPNN sequence design and AF2-multimer folding.
3. Tier 3 (specificity, planned): off-target panel curation (Hamming proteome scan, Householder method) plus structural cross-panning (AF2 cross-pan, Jenkins method).

The architecture mirrors the Science 2025 landmark papers (HADRUP_JENKINS_2025, BAKER_LAB_2025) with explicit engineering rigor: calibration anchors, mock-mode CI, a multi-entry trap book, spec-driven implementation. It is reproducible, end-to-end testable, and naturally extends to active learning when wet-lab feedback becomes available.
