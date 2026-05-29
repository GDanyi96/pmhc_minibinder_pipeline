# De novo pMHC-I minibinder design pipeline: scientific narrative

---

## Executive summary

This is a calibrated, reproducible funnel for designing small (~70–110 aa) de novo proteins that recognize a peptide-MHC class I complex — and, more consequentially, a **peptide-resolved interface metric** that measures whether a design reads the disease-defining peptide or merely the conserved MHC framework around it.

Across two structurally distinct design campaigns, that metric showed three things. First, the standard AF2 interface-confidence funnel is blind to the peptide-versus-framework distinction: designs that score in the positive-control band can be making zero contact with the peptide. Second, applied retroactively, the metric reclassified this project's own cycle-02 "hero" (design_2079_seq00) as a peptide-blind framework binder — its closest residue sits 28–40 Å from the peptide. Third, across 288 cycle-03 designs it surfaced a single design (design_3010_seq00) that contacts the WT1 peptide's specificity-bearing residues (N5, Y8) in the control-grade band, with one binder residue (R55) reading both.

This is a negative-result-rich pilot, and it is reported as one. Its value is methodological: the calibrated funnel, the metric audit that found two incompatible iPAE definitions across cycles, the peptide-contact decomposition that brings specificity information forward into the structure-prediction stage, an honest diagnosis of the limiting failure mode (framework bias is structural, not compositional — a charge hypothesis was raised and falsified), and a literature-anchored plan to fix it at scale. The architecture mirrors the platforms from Johansen et al. [1] and Liu et al. [2]. **All results in this writeup are in silico predictions; no experimental validation has been performed.**

---

## 1. Clinical motivation

The immune system surveils cells by reading peptides displayed on MHC class I — the cell's intracellular "library card." Tumor cells display peptides derived from oncogenic drivers and neoantigens this way. The natural recognition machinery (T-cell receptors) is hard to isolate, often low-affinity, and difficult to manufacture at scale.

De novo minibinders for pMHC-I complexes address several problems at once: they access intracellular antigens via the surface-presented peptide, which conventional antibodies cannot target; they are small, stable, and cheap to produce; and they can be built into chimeric antigen receptor (CAR) constructs to redirect T-cell killing — the "BIKE" concept set out in Johansen et al. [1]. A single computational pipeline can, in principle, address arbitrary tumor antigens. The broader clinical case across CAR-T, T-cell engagers, and autoimmune indications is laid out in Hadrup [3].

The bottleneck is **specificity**: the binder must recognize the cancer peptide but not closely related self-peptides presented on the same HLA on healthy cells. Because the MHC framework is conserved and the peptide occupies a small, partly buried patch, the dominant failure mode is a binder that engages the framework and ignores the peptide. That is the hard problem this pipeline is built to make measurable.

---

## 2. Target choice

The cognate complex is **WT1 RMFPNAPYL / HLA-A\*02:01**, crystal structure PDB 3HPJ.

| Component | Value | Rationale |
|---|---|---|
| Antigen | WT1 (Wilms Tumor 1) | Overexpressed in AML and multiple solid tumors; well-established TAA |
| Peptide | RMFPNAPYL (residues 126–134) | Top-tier cancer immunotherapy target |
| HLA allele | HLA-A\*02:01 | Most common allele in the Caucasian population (~50%); broadest applicability |
| Structure | PDB 3HPJ | Published crystal complex; clean starting point |

Liu et al. [2] targeted this exact complex, giving the project a direct experimental literature anchor: their published binders become positive controls; their results become the benchmark.

---

## 3. Pipeline architecture

Three sequential tiers — calibration, sensitivity, specificity — implemented as a Snakemake DAG with per-stage halt gates and mock-mode CI.

```
TIER 1  CALIBRATION   published controls through the assay, matched to target geometry
                      Q: do the metrics distinguish binders from non-binders?  A: yes
        ─────────────────────────────────────────────────────────────────────────────
TIER 2  SENSITIVITY   Stage 0 target prep → Stage 1 RFdiffusion backbones
                      → Stage 2 ProteinMPNN sequences → AF2-multimer fold + interface metrics
                      Cycle 02 (de novo) and Cycle 03 (BAKER-scaffold partial diffusion)
                      + the peptide-resolved iface_pep decomposition (this work)
        ─────────────────────────────────────────────────────────────────────────────
TIER 3  SPECIFICITY   AF2 cross-pan of survivors × off-target panel + Hamming proteome scan
                      (planned, cycle 04+)
```

The contribution of this work sits inside Tier 2: a decomposition of the AF2 interface metric into peptide and framework components, which pulls a specificity signal forward from Tier 3 into the per-design fold itself.

---

## 4. Tier 1: calibration via regime-matched controls

Before generating designs, the pipeline is run on published reference binders. AF2 confidence is not absolute — an iPAE of 8 Å is meaningless without empirical anchors. Two control sets were run, each **matched to the target geometry of the cycle it calibrates**:

- **Cycle-1 controls** (`results/cycle_01/stage2/metrics.json`) — full four-chain target (HLA heavy + β2m + peptide + binder), with TCR-sized published controls. Calibrates the cycle-02 de novo pipeline.
- **Cycle-3 truncated-baseline controls** — three-chain truncated target (HLA α1/α2 groove + peptide + binder), matching the BAKER-scaffold cycle-03 geometry.

| Control | Identity | What it rules out |
|---|---|---|
| P1 | Baker lab WT1 binder [2], cognate target | Thresholds too tight (a real positive must pass) |
| P2 | Jenkins NY1-B04 binder [1], cognate NY-ESO-1 | Calibration is WT1-only |
| P3 | Baker WT1 binder on the *wrong* target (MART-1) | Pipeline rewards generic foldedness |
| N1 | Scrambled Baker WT1 sequence, same backbone | Metrics are backbone-dominated / sequence-blind |
| N2 | Random 85-aa sequence | AF2 reports false-positive interfaces for anything |

**Cycle-1 (full target).** Positives P1/P2/P3 at iPAE 4.54–4.88 Å, ipLDDT 94–96; negatives N1/N2 at 24.7–25.7 Å, ipLDDT 31–35. Clean ~20 Å iPAE and ~60-point ipLDDT separation; ipLDDT is the sharper signal.

**Cycle-3 (truncated target), interface-8 Å decomposition.** Positives `iface_tot` 1.27–1.71 / `iface_pep` 1.46–2.60 / ipLDDT 93–97; negatives `iface_tot` 12.7–13.4 / `iface_pep` inf / ipLDDT 58–63.

Two points matter for everything downstream:

1. **The two control sets are not numerically comparable** (4.5–4.9 full-target vs 1.27–1.71 truncated). They calibrate different geometries. Cross-cycle magnitude comparisons of iPAE between a full-target design (cycle 02) and a truncated-target design (cycle 03) are therefore confounded by the target change, independent of any design improvement.
2. **The P3 result is the empirical justification for Tier 3.** P3 (Baker WT1 binder on MART-1) is a deliberate specificity failure. On overall metrics it scores like a true positive — the documented AF2 specificity blind spot, quantified by Mares et al. [5] at AUROC 0.06–0.22 for sequence predictors on structurally valid peptides. But in the interface-8 Å *decomposition*, P3's `iface_pep` degrades (2.60 vs P1's 1.70) while its `iface_mhc` stays intact (1.42). The peptide channel carries specificity information the aggregate metric hides. This observation motivates the metric used throughout (§7).

---

## 5. The metric audit

The two cycles stored their interface metrics under **different operational definitions of iPAE**, a discrepancy found and resolved during analysis (validated dataset-wide; recomputation matched stored values to Δ = 0.0000):

| | Cycle 02 (`metrics.jsonl`) | Cycle 03 (`manual_metrics`) |
|---|---|---|
| iPAE | **interface-8 Å**: mean `min(PAE_ij, PAE_ji)` over binder↔target pairs within 8 Å heavy-atom distance (Johansen et al. [1], Fig. 1B) | **position-slice**: mean PAE over *all* binder↔target pairs, no distance filter (Bennett-style `pae_interaction`) |
| ipLDDT | interface-contact-residue mean | binder-chain mean |
| Comparable across cycles? | only **ipTM** is directly comparable as stored | |

Both definitions are published and legitimate. Interface-8 Å is the better ranking metric for this system — it is size-robust (binder lengths span 54–118 aa) and is what the Hadrup/Jenkins platform reports [1]; the position-slice is the Bennett field default and is retained for baseline comparability. All candidate rankings in this writeup use the recomputed interface-8 Å definition consistently across both cycles. The slice/interface mismatch is also why the originally slice-ranked cycle-03 top list (§7) reordered substantially once recomputed.

---

## 6. Tier 2, Cycle 02: de novo generation and the reclassified hero

Cycle 02 ran the full four-chain target end to end with de novo RFdiffusion backbones, no ProteinMPNN amino-acid bias, and `num_recycles=3`. Full detail in [cycle_02.md](cycle_02.md).

| Stage | Input | Output | Yield |
|---|---|---|---|
| RFdiffusion backbones | — | 200 | — |
| Geometry pass | 200 | 26 | 13% |
| ProteinMPNN (4 seqs/backbone) | 26 | 104 | — |
| AF2-multimer fold (top 50) | 50 | 50 | — |
| Strict cut (iPAE ≤ 12 **and** ipLDDT ≥ 88) | 50 | **1** | **2%** |

The interface-8 Å distribution is informative: median 22.66 Å, ~70% of designs (35/50) in the 20–25 Å range (the cycle-1 negative zone), min 6.41 Å (the single passing design). AF2 reads most designs as coherent four-helix bundles that simply do not contact the binding cleft with confidence — the placement deficit propagating from Stage 1, which ProteinMPNN cannot repair.

**The hero, reclassified.** design_2079_seq00 (99 aa, four-helix bundle, BSA 976 Å² in the Liu et al. [2] minibinder range) crossed the strict cuts at iPAE 6.41 Å / ipLDDT 91.1, and the original narrative reported it as a positive-band win. The peptide-contact decomposition (§7) overturns that reading:

| Metric | Value | Reading |
|---|---|---|
| `iface_tot` | 6.41 | Above the cycle-1 full-target positive band (4.5–4.9) — a weak interface even in aggregate |
| `iface_pep` | **inf** | **No binder atom within 8 Å of any peptide atom** |
| `iface_mhc` | 6.41 | The entire interface is MHC framework |
| closest binder→peptide distance | **28–40 Å** | The binder docks on a distal MHC surface, nowhere near the groove |

design_2079 is a peptide-blind MHC-framework binder. Its 6.41 was entirely framework contact; the aggregate metric never revealed that the peptide — the entire point of the target — was 28–40 Å away. Two further qualifications from the original analysis stand: the backbone is a sequence-design lottery (its four ProteinMPNN sequences scored 6.41 / 19.26 / 24.11 / 24.93 — only one finds a binder-like configuration), and the sequence is 40.4% alanine (the unconstrained ProteinMPNN regularizer trap; AF2-orthogonal but wet-lab-consequential, with no W/Y for UV quantitation). The alanine issue was fixed in cycle 03; the peptide-blindness was not, because it is upstream of sequence design.

---

## 7. The peptide-resolved interface metric

The standard de novo binder funnel ranks on ipTM and an aggregate interface iPAE. For a generic protein target that is sufficient. For pMHC it is not: the target surface is mostly conserved framework with a small variable peptide patch, and a binder can satisfy every aggregate metric while contacting only the framework — exactly the cycle-02 hero, and (§8) most of cycle 03.

The fix is to decompose the interface-8 Å iPAE by target sub-chain:

- `iface_pep` — restricted to binder↔peptide residue pairs within 8 Å
- `iface_mhc` — restricted to binder↔(HLA + β2m) pairs within 8 Å
- `iface_tot` — all binder↔target pairs

`iface_pep = inf` is a hard, interpretable signal: the binder makes no contact with the peptide at all. The control panel calibrates the metric directly — positives sit at `iface_pep` 1.46–2.60, and the P3 specificity-failure control degrades on this channel while holding on `iface_mhc` (§4). This is the methodological core of the project: a specificity signal computed from a single complex prediction, before any cross-pan, that standard funnels discard by aggregation. It is an in-silico signal, not an affinity measurement.

---

## 8. Cycle 03 (scaffold-based) — summary

Cycle 03 (sub-run A) replaced de novo generation with partial diffusion (`partial_T=15`) from 152 aligned BAKER scaffolds onto the truncated three-chain target, with ProteinMPNN amino-acid bias (A: −2.0; E/L/R: +1.0) and `num_recycles=6`. Full detail in [cycle_03.md](cycle_03.md); the headline results:

- **Funnel:** 150 backbones (150 unique scaffolds) → 72 geometry-pass (48%, up from 13%) → 288 sequences → 288 AF2 folds.
- **Peptide engagement:** only 91/288 (32%) make any peptide contact; median engager `iface_pep` 14.2 (barely-contact); **one** design in the positive-control band.
- **Top by standard metrics:** design_3084_seq02 (`iface_tot` 1.56, ipTM 0.89 — the best overall interface, exactly what a standard funnel would rank first) and design_3054_seq00 (`iface_tot` 1.61, `iface_pep` inf — confirmed framework-only). Neither is the peptide reader; a standard funnel picks one of them and never sees 3010.
- **The one reader:** design_3010_seq00 (§9).
- **Scaffold lineage finding:** cross-allele scaffolds (HLA-CA RMSD 4–6 Å) passed at 54% vs 38% for A\*02:01-native (<1 Å) scaffolds — HLA structural similarity is *not* predictive of transfer quality (see [methodological_lessons.md](methodological_lessons.md)).

---

## 9. The validated lead: design_3010_seq00

One design in 288 contacts the peptide in the control-grade band. From scaffold scaf109 (cross-allele, BAKER cluster 6):

| Metric | design_3010_seq00 | Positive band | Framework champ (3054) |
|---|---|---|---|
| `iface_tot` | 2.98 | 1.27–1.71 | 1.61 |
| `iface_pep` | **2.10** | 1.46–2.60 | inf |
| `iface_mhc` | 3.35 | — | — |

3010 is the only design that reads the peptide more tightly than the framework. Per-residue contact analysis (AF2-predicted geometry; chains confirmed A=HLA / B=peptide RMFPNAPYL / C=binder, no numbering offset):

- Tight (≤ 5 Å) contact across the entire solvent-exposed central bulge **P4–P8**, including **both** specificity residues — N5 at 2.34 Å and Y8 at 2.71 Å (Y8 PAE 4.17, a confident contact).
- The A\*02:01 anchors (P2-Met, P9-Leu) are not tightly engaged (P2 at 6.81 Å; P9 grazes at 4.40 Å) — the binder reads the variable, specificity-bearing face, not the buried anchors.
- A contiguous binder segment (residues 48–60) forms the reader. **R55 contacts both N5 and Y8** — the specificity linchpin; E57 reads P4/P5.

This is the "reads WT1" outcome rather than the "grazes the flanks" outcome, and it yields a falsifiable, residue-level prediction: substituting N5 or Y8, or mutating binder R55, should collapse `iface_pep`. That is the decisive cycle-04 specificity test (§11). **Caveats:** n = 1; in silico only; and 3010 is a weaker overall binder than the validated positives (`iface_tot` 2.98 vs 1.27–1.71). It is a characterized lead, not a result.

---

## 10. The central finding: framework bias is structural

The limiting failure mode across both cycles is that designs engage the conserved MHC framework rather than the variable peptide. This is the field's central pMHC challenge [1][2][4], and the pipeline made it quantitative: cycle-02 peptide engagement 40%, cycle-03 32%, with one control-grade peptide reader total.

A compositional explanation was hypothesized and tested. The cycle-03 amino-acid bias drove net charge strongly anionic (−1.6 → −9.1), and the MHC α1/α2 walls present a large charged surface while the peptide is small and partly buried — suggesting charged binders find electrostatic complementarity with the framework. **The hypothesis was falsified.** Peptide-engaging designs were *more* anionic than peptide-blind ones (mean net charge −12.4 vs −7.6, the opposite of the predicted direction), with a weak reverse within-engager correlation (~0.11). The decisive evidence is single-backbone: on scaffold scaf158, seq00 gives `iface_pep` inf while seq01 gives 5.23 — same backbone, different sequence, peptide contact modulated only at the margin. Backbone geometry sets the regime; sequence cannot rescue a backbone whose binder is placed away from the peptide. Peptide-blindness is structural — an RFdiffusion placement property — not compositional. The fix therefore belongs at the generation stage, not in sequence design.

Including MHC-framework hotspots was not an error: the validated positive controls P1/P2 contact *both* peptide and framework (a 9-mer alone is too small an interface for affinity). The empirical finding is that the *balance* skewed entirely to framework — and that BAKER-scaffold reuse, by inheriting framework-favoring topologies, slightly worsened it relative to de novo generation.

---

## 11. Honest cross-cycle comparison

Cycles 02 and 03 are **not** a controlled A/B test. At least five knobs changed simultaneously: generation mode (de novo → partial diffusion from BAKER scaffolds), ProteinMPNN bias (none → A/E/L/R bias), AF2 recycles (3 → 6), target geometry (full four-chain → truncated three-chain), and the stored iPAE definition (interface-8 Å → position-slice). No single-variable effect can be attributed.

What survives the confound is the qualitative, geometry-independent finding: **both** strategies are dominated by peptide-blind framework binders (cycle-02 hero `iface_pep` inf; cycle-03 framework champions `iface_pep` inf; 32–40% peptide engagement; one reader). The interface-8 Å recomputation removes the metric-definition confound but not the target-geometry one — so the cleanest cross-cycle statement is that framework bias is robust to the generation strategy, not that one cycle's binders are "N× better" than the other's. A single forced synthesis claiming one effect would misrepresent a five-knob change; this comparison is reported as the bounded statement the data support.

---

## 12. Cycle 04 plan

The diagnosis (placement, not composition) and the current tool landscape define the next campaign. The generative model is not the lever — running more RFdiffusion with the same conditioning would produce more framework binders.

1. **Fix conditioning (the lever).** Peptide-only hotspot conditioning at Stage 1 — drop the MHC-framework hotspots (A65/66/150/155), keep peptide hotspots — so the diffusion trajectory anchors on the peptide. Free; directly targets the diagnosed cause.
2. **Scale 10–50×.** The field tests thousands of designs to validate a handful; this pilot tested one. Scale the front of the funnel to match.
3. **Generate in parallel with two engines.** RFdiffusion (for topological diversity and the partial-diffusion control already working) **and** BindCraft [7], an AF2-hallucination pipeline reporting 10–100% experimental success (avg ~46%) with fewer than 100 designs needing testing, outperforming RFdiffusion, AlphaProteo, and Chai-2. BindCraft is not a pMHC solution by itself — it targets a surface, so the peptide-versus-framework competition persists, and it has reported failures (PD-1/PD-L1) — but it is the strongest current option for raw hit rate. Note that RFdiffusion2 [9] is an enzyme/atomic-motif scaffolding model, *not* a binder upgrade; the binder-relevant RFdiffusion2 work is a separately fine-tuned antibody network [10] (a different modality, though notably applied to a Phox2b pMHC scFv).
4. **Upgrade the filter cascade.** A meta-analysis of 3,766 experimentally tested binders found interface-focused metrics — notably the AF3-derived ipSAE — outperform ipTM and `pae_interaction` for predicting experimental success [8]. Cascade: ipSAE + AF3/Boltz-2 re-prediction + the `iface_pep` peptide-contact gate (this work) + short MD on survivors (Johansen et al. [1] used MD to catch failures that passed ipTM/pLDDT).
5. **Decisive specificity test for 3010.** Fold 3010 against single-position peptide mutants (Y8A, N5A) and the IEDB off-target panel; if 3010 genuinely reads N5/Y8 via R55, `iface_pep` should collapse while `iface_mhc` holds — the in-silico analog of the Bentzen et al. [6] fingerprinting alanine scan.
6. **Scaffold pre-filtering (from cycle-03 findings).** Pre-select scaffolds by binder-to-peptide proximity rather than HLA-CA RMSD, by binder length (to remove the length-out-of-range failures), and from high-pass clusters.

---

## 13. Tier 3 and beyond (planned)

**Tier 3 — specificity.** Sequence-based off-target panel curation by Hamming distance to RMFPNAPYL [4] (CPU-only, peptide-dependent only), then AF2/Boltz cross-panning of survivors against each off-target pMHC, using `iface_pep` as the contrast axis. The Jenkins platform [1] performs the equivalent cross-pan before any wet-lab handoff.

**Stage 4 — diversity.** ESM-2 embedding + farthest-point sampling to select a structurally diverse subset for synthesis.

**Stage 5 — active learning.** A surrogate (LightGBM / GP, UCB acquisition) trained on design–outcome pairs once wet-lab feedback exists — the experimental partnership a PhD project would build, with Bentzen-style DNA-barcoded MHC-multimer fingerprinting [6] as the canonical readout.

---

## 14. Reproducibility

The repo bootstraps on a RunPod A100 pod via `bootstrap.sh` (two Python environments, pinned JAX/PyTorch, model weights). Snakemake mock mode runs the full DAG on synthetic fixtures in under a second and gates CI on every push. All numerical thresholds are externalized to `configs/thresholds.yaml`; no magic numbers. Per-stage halt gates surface upstream breakage as loud failures rather than silent fallbacks. Engineering-level gotchas are logged separately in `docs/known_traps.md`; the conceptual pitfalls that bear on the science are in [methodological_lessons.md](methodological_lessons.md). Heavy run outputs (raw AF2 predictions, the full design set) stay on the pod and are reproducible from the committed configs; only the curated analysis tables, the one validated structure (design_3010), and the two hero contact tables are committed.

---

## 15. References

[1] Johansen, K.H., Wolff, D.S., Scapolo, B., Fernández-Quintero, M.L., Christensen, C.R., Loeffler, J.R., Rivera-de-Torre, E., Overath, M.D., Munk, K.K., Morell, O., Viuff, M.C., Lacunza, I., Damm Englund, A.T., Due, M., Gharpure, A., Forli, S., Rodriguez Pardo, C., Tamhane, T., Andersen, E.Q., Björnsson, K.H., Fernandes, J.S., Voss, L.F., Thumtecho, S., Ward, A.B., Ormhøj, M., Hadrup, S.R., Jenkins, T.P. De novo-designed pMHC binders facilitate T cell-mediated cytotoxicity toward cancer cells. *Science* **389**(6758) (2025). DOI: [10.1126/science.adv0422](https://doi.org/10.1126/science.adv0422). PDB: 9NNF.

[2] Liu, B., Greenwood, N.F., Bonzanini, J.E., Motmaen, A., Meyerberg, J., Dao, T., Xiang, X., Ault, R., Sharp, J., Wang, C., Visani, G.M., Vafeados, D.K., Roullier, N., Nourmohammad, A., Scheinberg, D.A., Garcia, K.C., Baker, D. Design of high-specificity binders for peptide-MHC-I complexes. *Science* **389**, 386 (2025). DOI: [10.1126/science.adv0185](https://doi.org/10.1126/science.adv0185). Code: [Zenodo (77forest/pmhci_binder_design)](https://doi.org/10.5281/zenodo.15169815).

[3] Hadrup, S.R. Artificial intelligence is expediting the development of therapeutic immunotherapies. *ESMO Daily Reporter*, ESMO Immuno-Oncology Congress 2025 opinion (2025). Keynote: "Using AI to advance therapeutic development of immunotherapies." Available at [dailyreporter.esmo.org](https://dailyreporter.esmo.org/esmo-immuno-oncology-congress-2025/opinions/artificial-intelligence-is-expediting-the-development-of-therapeutic-immunotherapies).

[4] Householder, K.D., Xiang, X., Jude, K.M., Deng, A., Obenaus, M., Zhao, Y., Wilson, S.C., Chen, X., Wang, N., Garcia, K.C. De novo design and structure of a peptide-centric TCR mimic binding module. *Science* **389**(6758), 375–379 (2025). DOI: [10.1126/science.adv3813](https://doi.org/10.1126/science.adv3813). PDB: 9MIN.

[5] Mares, S.E., Espinoza Weinberger, A., Ioannidis, N.M. Generation of structure-guided pMHC-I libraries using Diffusion Models. *2nd International Conference on Machine Learning in Generative AI and Biology Workshop* (2025). Code: [github.com/sermare/struct-mhc-dev](https://github.com/sermare/struct-mhc-dev).

[6] Bentzen, A.K., Such, L., Jensen, K.K., Marquard, A.M., Jessen, L.E., Miller, N.J., Church, C.D., Lyngaa, R., Koelle, D.M., Becker, J.C., Linnemann, C., Schumacher, T.N.M., Marcatili, P., Nghiem, P., Nielsen, M., Hadrup, S.R. T cell receptor fingerprinting enables in-depth characterization of the interactions governing recognition of peptide-MHC complexes. *Nature Biotechnology* **36**(12), 1191–1196 (2018). DOI: [10.1038/nbt.4303](https://doi.org/10.1038/nbt.4303).

[7] Pacesa, M., Nickel, L., Schmidt, J., Pyatova, M., Schellhaas, C., Kissling, L., Sankaran, S., Ahmed, T., Bonati, J., Rosset, S., Wang, C., Dauparas, J., Ovchinnikov, S., Correia, B.E. One-shot design of functional protein binders with BindCraft. *Nature* (2025). DOI: [10.1038/s41586-025-09429-6](https://doi.org/10.1038/s41586-025-09429-6).

[8] Overath, M.D., Rygaard, K.B., et al. Predicting experimental success in de novo binder design: a meta-analysis of 3,766 experimentally characterised binders. Preprint (2025). Introduces the AF3-derived ipSAE as a target-agnostic interface filter outperforming ipTM and `pae_interaction`.

[9] Ahern, W., Yim, J., Tischer, D., Salike, S., Woodbury, S., Kim, D., Kalvet, I., Kipnis, Y., Coventry, B., Altae-Tran, H., Bauer, M., Barzilay, R., Jaakkola, T., Krishna, R., Baker, D. Atom-level enzyme active site scaffolding using RFdiffusion2. *Nature Methods* (2025). DOI: [10.1038/s41592-025-02975-x](https://doi.org/10.1038/s41592-025-02975-x).

[10] Bennett, N.R., et al. Atomically accurate de novo design of antibodies with RFdiffusion. *Nature* (2025). DOI: [10.1038/s41586-025-09721-5](https://doi.org/10.1038/s41586-025-09721-5).

---
**Committed data and code:** the cross-cycle funnel
[results/master_design_journey.csv](../results/master_design_journey.csv); the cycle-03
artifacts in [results/cycle_03/analysis/](../results/cycle_03/analysis/) (candidate
dossier, scaffold lineage, the validated design_3010 structure); and the per-residue
contact analyses [12_design3010_peptide_contacts.py](../analysis/scripts/12_design3010_peptide_contacts.py)
and [13_design2079_peptide_contacts.py](../analysis/scripts/13_design2079_peptide_contacts.py).


## TL;DR for the application

A calibrated, reproducible de novo pMHC-I minibinder pipeline mirroring the Johansen et al. [1] and Liu et al. [2] *Science* 2025 platforms, distinguished by a **peptide-resolved interface metric** that measures whether a design reads the disease-defining peptide or the conserved MHC framework. That metric (a) exposed a blind spot in the standard AF2 funnel, (b) found and reconciled two incompatible iPAE definitions across cycles, (c) reclassified this project's own cycle-02 hero as a peptide-blind framework binder (closest residue 28–40 Å from the peptide), (d) surfaced one design (3010) that reads both WT1 specificity residues N5/Y8 via binder R55 in the control-grade band, and (e) supported an honest, falsified charge hypothesis that localized the dominant failure (framework bias) to RFdiffusion placement rather than sequence composition. The cycle-04 plan fixes the diagnosed cause with peptide-only conditioning, scales 10–50×, runs RFdiffusion and BindCraft [7] in parallel, and adds an ipSAE-based filter cascade [8]. All results are in silico; no experimental validation has been performed.
