# Methodological lessons

Conceptual pitfalls in de novo pMHC-I binder design surfaced by this project, each with the evidence and the design consequence. These are the science-bearing lessons; purely engineering gotchas (environment, dependency, and I/O issues) are logged separately in `known_traps.md` and kept out of the scientific narrative. All evidence is in-silico AF2-predicted unless noted.

---

## 1. AF2 interface confidence is blind to peptide-versus-framework engagement

A design can score in the positive-control band on aggregate interface iPAE and ipTM while making **zero** contact with the peptide. The cycle-02 hero (design_2079_seq00) passed the strict cuts at iPAE 6.41 / ipLDDT 91 — and its closest residue is 28–40 Å from the peptide; the entire interface is MHC framework.

- **Evidence:** decomposing interface-8 Å iPAE into `iface_pep` / `iface_mhc` gives `iface_pep = inf` (no binder atom within 8 Å of the peptide) for the cycle-02 hero and for the cycle-03 framework champions, despite passing every aggregate metric.
- **Consequence:** rank pMHC designs on a peptide-restricted interface metric, not on aggregate ipTM/iPAE. The peptide channel is the specificity axis; the aggregate hides it.

---

## 2. "iPAE" names two non-interchangeable metrics

Cycle 02 and cycle 03 stored interface confidence under different operational definitions of the same field name. Comparing them as if identical mis-ranks candidates.

- **Definitions:** interface-8 Å iPAE = mean `min(PAE_ij, PAE_ji)` over binder↔target pairs within 8 Å (Johansen et al. [1], Fig. 1B); position-slice iPAE = mean PAE over all binder↔target pairs, no distance filter (Bennett-style `pae_interaction`).
- **Evidence:** of the slice-ranked cycle-03 top 10, only 3 engage the peptide; recomputing on interface-8 Å reordered the list substantially (recomputation matched stored values to Δ = 0.0000, confirming the difference is definitional, not a bug).
- **Consequence:** fix one ranking definition and apply it consistently across cycles. Interface-8 Å is size-robust (binder lengths 54–118 aa) and matches the Hadrup/Jenkins platform; the slice is the Bennett default and is kept only for baseline comparability.

---

## 3. Controls must be regime-matched to the target geometry

Calibration numbers do not transfer across target constructs. Full-target controls (cycle 1) sit at iPAE 4.5–4.9 (positives); truncated-target controls (cycle 3) sit at `iface_tot` 1.27–1.71 — the same published binders, different geometry.

- **Consequence:** run a control set per target geometry, and never compare a full-target design's iPAE magnitude directly to a truncated-target design's. The target change is a confound in any cross-cycle magnitude comparison (see narrative §11) — the only safe cross-cycle statement here is qualitative (both cycles framework-biased).

---

## 4. HLA structural similarity does not predict scaffold transfer quality

The intuitive prior — a scaffold whose HLA aligns tightly to the target transfers cleanly — is false and backwards.

- **Evidence:** A\*02:01-native scaffolds (HLA-CA RMSD < 1 Å) passed at 38%; cross-allele scaffolds (4–6 Å) passed at 54%. All five top candidates came from cross-allele scaffolds.
- **Mechanism:** native scaffolds inherit their original target-peptide register and `partial_T=15` cannot remodel them onto WT1; cross-allele scaffolds are less constrained and get repositioned toward the groove by alignment.
- **Consequence:** pre-filter a scaffold library on binder-to-peptide proximity, not HLA-HLA RMSD.

---

## 5. Peptide-blindness is structural, not compositional

A tempting compositional explanation — strongly anionic binders (cycle-03 net charge −9) favoring the charged MHC α1/α2 walls over the small, partly buried peptide — is wrong.

- **Evidence:** peptide-engaging designs were *more* anionic than peptide-blind ones (mean net charge −12.4 vs −7.6), the opposite of the prediction; within-engager correlation was weak and reversed (~0.11). The decisive control is single-backbone: on scaf158, seq00 gives `iface_pep` inf while seq01 gives 5.23 — same backbone, different sequence, peptide contact barely moves.
- **Consequence:** the fix for framework bias belongs at the generation stage (RFdiffusion placement / hotspot conditioning), not in sequence design. ProteinMPNN cannot rescue a backbone whose binder is placed away from the peptide.

---

## 6. Multi-knob cycles are not controlled comparisons

Cycles 02 and 03 changed five things at once (generation mode, MPNN bias, AF2 recycles, target geometry, iPAE definition). No single-variable effect is attributable.

- **Consequence:** report the bounded statement the data support — framework bias is robust across both de novo and scaffold-based generation — rather than a single-cause synthesis. A future cycle that isolates one knob (e.g. peptide-only vs mixed hotspots, all else fixed) is required for a causal claim.

---

## 7. Sequence predictors and aggregate structure scores both miss novel-peptide specificity

The deliberate specificity-failure control (P3: a WT1 binder folded against MART-1) scores like a true positive on aggregate metrics — the documented AF2 specificity blind spot, quantified by Mares et al. [5] at AUROC 0.06–0.22 for sequence predictors on structurally valid peptides.

- **Evidence in this pipeline:** P3's `iface_pep` degrades (2.60 vs P1's 1.70) while `iface_mhc` holds (1.42) — the decomposed metric recovers the specificity signal the aggregate loses.
- **Consequence:** treat sequence-based pMHC predictors (NetMHCpan/MHCflurry) as baselines, not oracles, and assess specificity structurally — via the peptide-restricted interface channel here, and via explicit cross-panning (Tier 3) before any wet-lab handoff.

---

*References as numbered in [narrative.md](narrative.md).*
