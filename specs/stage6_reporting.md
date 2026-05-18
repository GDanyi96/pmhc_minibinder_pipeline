# Stage 6 — Cycle Reporting

## Goal

Aggregate per-cycle artefacts into a single publication-grade markdown
report with figures, control-panel verification, and a decision log.
The report is committed to the repo (in `reports/cycle_NN.md`) so a
reviewer can read the project's progress at a glance.

## Inputs

- `results/cycle_NN/stage2/metrics.parquet` + `controls_metrics.json`.
- `results/cycle_NN/stage3/specificity_scores.parquet`.
- `results/cycle_NN/stage4/diversity_report.json`.
- `results/cycle_NN/stage5/al_metrics.json` +
  `next_cycle_candidates.json`.
- `results/cycle_NN/HALT.txt` (if present).

## Outputs

- `reports/cycle_NN.md` — narrative report with embedded figures.
- `reports/cycle_NN_assets/` — PNG/SVG figures, committed.
  - `fig01_pipeline_funnel.png` — N designs surviving each stage.
  - `fig02_ipae_distribution.png` — iPAE histogram with control overlays.
  - `fig03_crosspan_heatmap.png` — ΔiPAE per (design, off-target).
  - `fig04_diversity_tsne.png` — t-SNE of ESM embeddings, FPS-selected
    highlighted.
  - `fig05_al_calibration.png` — surrogate predicted vs actual; ranking.

## Tools

- `matplotlib`, `seaborn` for figures.
- `pandas` for table aggregation.
- `jinja2` for markdown templating.

## Anchored references

- `HADRUP_JENKINS_2025` Fig 1 — figure style reference for funnel +
  cross-panning heatmap.

## Implementation tasks

- [ ] `render_report.py`: load all stage outputs into a `CycleReport`
      pydantic model.
- [ ] Render funnel figure: bar chart of survivors after each stage.
- [ ] Render iPAE histogram with vertical lines for P1–N2 controls.
- [ ] Render cross-pan heatmap (designs × off-targets, colour = ΔiPAE).
- [ ] Render t-SNE of ESM embeddings, FPS picks highlighted.
- [ ] Render AL calibration plot + top-K acquisition table.
- [ ] Fill markdown template; commit figures + report.
- [ ] If `HALT.txt` exists, prepend a HALT section explaining which
      control failed and what to inspect.
- [ ] `--mock`: cp `tests/fixtures/report/sample_manifest.json` to output.

## Verification criteria

- Smoke test: `render_report.py --mock` exits 0; writes a stub
  manifest.
- Real test: cycle 1 report renders end-to-end with all 5 figures;
  passes `markdownlint`.
- Reviewer test: a non-domain-expert can read the report in 5 minutes
  and explain what the cycle achieved.

## Pitfalls

- Figures must be reproducible: seed numpy/matplotlib RNGs.
- Commit PNG **and** SVG for each figure (SVG for resolution-independent
  inclusion in a PDF; PNG for GitHub preview).
- The HALT path needs to be tested in mock mode at least once.
- Do not put PII / pod paths in the committed report.
