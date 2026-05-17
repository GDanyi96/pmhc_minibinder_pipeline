# Stage 3 — In Silico Cross-Reactivity Panning

## Goal

For every survivor from stage 2, re-fold the binder against an off-target
peptide panel presented on HLA-A\*02:01 and compute a specificity score
ΔiPAE = iPAE(off-target) − iPAE(on-target). Reject designs that bind
off-target peptides almost as well as the on-target. This is the
specificity bottleneck — the differentiating depth of this pipeline.

## Inputs

- `results/cycle_NN/stage2/survivors.json` — designs passing stage 2.
- `results/cycle_NN/stage2/colabfold/design_*/ranked_0.pdb` — on-target structures.
- `configs/thresholds.yaml` (`crosspan.panel`) — off-target peptide list:
  MART-1, gp100, NY-ESO-1, MAGE-A3, CMV-pp65, Flu-M1.
- `data/targets/3hpj_clean.pdb` — MHC scaffold for peptide swapping.

## Outputs

- `results/cycle_NN/stage3/crosspan_matrix.parquet` — rows: designs,
  columns: off-target peptides, values: iPAE.
- `results/cycle_NN/stage3/specificity_scores.parquet` — per-design min
  ΔiPAE across panel.
- `results/cycle_NN/stage3/specificity_survivors.json` — designs passing
  `crosspan.min_delta_ipae`.

## Tools

- Same ColabFold image as stage 2; same `num_recycles=6`.
- `pandas` for matrix aggregation.

## Anchored references

- `HADRUP_JENKINS_2025` — Fig 1B cross-panning grid (real experimental
  validation; this stage is the in silico analog).
- `HOUSEHOLDER_GARCIA_2025` — proteome-wide Hamming distance scan as an
  upstream cheap filter (optional v2 enhancement; not in cycle 1).
- `BENTZEN_HADRUP_2019` — context: DNA-barcoded multimer experimental
  readout this stage mimics in silico.

## Implementation tasks

- [ ] `crosspan.py`: for each (survivor design, off-target peptide), build
      a hybrid PDB: clean MHC scaffold + swapped peptide + the design's
      predicted binder coordinates as starting structure.
- [ ] Re-fold the hybrid with AF2-multimer (same chain order, same recycles).
- [ ] Parse iPAE per (design, off-target). Assemble crosspan_matrix.
- [ ] Compute ΔiPAE per (design, off-target); take min across panel as the
      specificity score.
- [ ] Filter: reject designs where min ΔiPAE < `crosspan.min_delta_ipae` (3 Å).
- [ ] Run P3 control (Baker WT1 binder vs MART-1) and check it lands as
      expected (iPAE > P1 + 3 Å).
- [ ] `--mock`: cp `tests/fixtures/crosspan/sample_panel.json` to output.

## Verification criteria

- Smoke test: `crosspan.py --mock` exits 0 in <1s.
- Real test (pod): 6-peptide panel × 96 survivors = 576 AF2 jobs; should
  complete in <8 h on A100 with ColabFold remote MSA cache.
- P3 control passes (iPAE separation > 3 Å vs P1 on-target).

## Pitfalls

- Peptide swapping must preserve MHC anchor residue interactions; the
  scaffold (HC, β2m) is held fixed and peptide is replaced atom-by-atom
  using a backbone superposition. Side chains rebuilt by AF2.
- Off-target peptides must be the same length as the on-target (or use a
  separate length-stratified panel). WT1 is a 9mer; all panel peptides are
  9mers except MART-1 (10mer) — handle MART-1 with a separate
  HLA-A\*02:01 9mer alternative or use the canonical 10mer with care.
- Remote MSA cache: re-folding identical MHCs reuses cached MSAs and
  speeds up by ~5×. Don't disable caching.
- iPAE definition must match stage 2 exactly (binder ↔ {HC, β2m, peptide}).
