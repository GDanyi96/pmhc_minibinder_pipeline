# Stage 5 — Active Learning

## Goal

Train a surrogate model on (initially synthetic) experimental labels and
use it to score candidate designs for the next cycle, ranking by a
UCB acquisition function. LightGBM is the primary surrogate; a GP is
trained in parallel as a principled-uncertainty baseline. In cycle 1
labels are synthesised from the in silico metrics; cycle ≥ 2 ingests
real wet-lab labels (when available).

## Inputs

- `results/cycle_NN/stage4/fps_selection.json` + `embeddings.npz` —
  candidates for next cycle.
- `results/cycle_NN/stage2/metrics.parquet` — in silico labels for cold-start.
- `results/labels/cycle_NN.csv` (cycle ≥ 2) — real experimental labels.

## Outputs

- `results/cycle_NN/stage5/surrogate_lgbm.txt` — trained LightGBM model.
- `results/cycle_NN/stage5/surrogate_gp.pkl` — trained GP baseline.
- `results/cycle_NN/stage5/predictions.parquet` — mean + variance per
  candidate.
- `results/cycle_NN/stage5/al_metrics.json` — train/val MAE/R², calibration
  metrics, top-K acquisition rankings.
- `results/cycle_NN/stage5/next_cycle_candidates.json` — top-K by UCB.

## Tools

- `lightgbm` for primary surrogate.
- `scikit-learn` (`GaussianProcessRegressor`) for baseline.
- `numpy` for UCB.
- ESM-2 embeddings from stage 4 as input features (concat with engineered
  features: length, charge, hydrophobicity, iPAE, ipLDDT).

## Anchored references

- General active learning literature; no specific paper anchor for the
  acquisition strategy. Document the choice (UCB κ=2.0) in code.
- `HADRUP_JENKINS_2025` — labels at cycle ≥ 2 are barcoded multimer
  enrichment scores; cycle 1 uses synthetic labels derived from iPAE.

## Implementation tasks

- [ ] `active_learning.py`: assemble feature matrix (ESM embeddings +
      engineered features) from inputs.
- [ ] Synthesise cycle-1 labels: label = a · (-iPAE) + b · ipLDDT + noise.
      Document the synthetic-label function explicitly in the docstring.
- [ ] Train LightGBM with 5-fold CV; log MAE, R², calibration.
- [ ] Train GP baseline (kernel=Matérn 5/2) on the same split.
- [ ] Predict mean + variance for all FPS candidates.
- [ ] UCB acquisition: `score = μ + κ · σ`, κ from thresholds.yaml.
- [ ] Emit al_metrics.json + next_cycle_candidates.json.
- [ ] `--mock`: cp `tests/fixtures/active_learning/sample_predictions.json`.

## Verification criteria

- Smoke test: `active_learning.py --mock` exits 0 in <1s.
- Real test: cycle 1 trains on 96 synthetic labels in <30 s on CPU;
  predicted ranking has Spearman > 0.5 vs in silico iPAE (since labels are
  derived from iPAE).
- Cycle ≥ 2: train MAE decreases vs cycle 1 (active-learning gain).

## Pitfalls

- Cycle 1 synthetic labels are a calibration step, not a real result.
  Make this **loud** in the docstring and the cycle 1 report.
- ESM embeddings are 1280-D; engineered features must be normalised or
  LightGBM will overweight them.
- GP scaling: GaussianProcessRegressor is O(n³); cap training set at 500
  for the baseline.
- UCB κ=2.0 is reasonable for ~96 candidates; for very small batches
  (<10) prefer expected-improvement.
- Don't leak test-set designs into next_cycle_candidates.
