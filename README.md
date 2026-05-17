# pMHC-I Minibinder Design Pipeline

> De novo minibinder design against peptide-MHC class I targets, with built-in
> in silico cross-reactivity panning and an active-learning loop.
> Built as a portfolio project for the Jenkins / Hadrup lab at DTU.

**Primary target:** WT1 RMFPNAPYL / HLA-A\*02:01 (PDB 3HPJ) — a leukemia-relevant
tumour antigen.

**Why it matters:** The current bottleneck in TCR-mimic and minibinder design
is *specificity*, not affinity. This pipeline operationalises an in silico
cross-reactivity panel (stage 3) and an ESM-2-driven diversity gate
(stage 4) on top of an otherwise standard RFdiffusion → ProteinMPNN → AF2
stack.

## Pipeline

```
Stage 0  Target prep        — clean PDB, define hotspots
Stage 1  RFdiffusion        — generate de novo binder backbones
Stage 2  ProteinMPNN + AF2  — design sequences, validate folds (iPAE/ipLDDT)
Stage 3  Cross-pan          — score every survivor vs off-target peptide grid
Stage 4  Diversity          — ESM-2 embeddings + farthest-point sampling
Stage 5  Active learning    — LightGBM + GP surrogate; UCB acquisition
Stage 6  Reporting          — publication-grade cycle report
```

All stages are orchestrated by Snakemake; thresholds live in
`configs/thresholds.yaml`; per-stage contracts live in `specs/`.

## Quick start (on a GPU pod with `/workspace/pipeline/` checkout)

```bash
bash bootstrap.sh                              # installs deps, pulls Docker images, downloads weights
snakemake --dry-run --config mock=true -j1     # validates the DAG (no GPU needed)
snakemake --cores all                          # real cycle-1 run (A100 recommended)
```

## Repo layout

See [`INDEX.md`](INDEX.md) for a one-line description of every tracked file.
See [`CLAUDE.md`](CLAUDE.md) for the project rulebook and locked decisions.

## Scientific anchors

| Anchor                    | Role                                                      |
|---------------------------|-----------------------------------------------------------|
| `HADRUP_JENKINS_2025`     | Primary pipeline reference (RFdiffusion config, AF2 cuts) |
| `BAKER_LAB_2025`          | WT1 binder sequence, partial-diffusion scaffold reuse     |
| `HOUSEHOLDER_GARCIA_2025` | NY-ESO-1 alt binder, proteome Hamming-distance safety     |
| `MARES_IOANNIDIS_2025`    | Predictor blind spots; AF2 num_recycles=6                 |
| `BENTZEN_HADRUP_2019`     | Wet-lab cross-reactivity readout (motivation only)        |

## Status

Cycle 1 in progress. Per-cycle reports land in `reports/cycle_NN.md`.

## License

MIT — see [`LICENSE`](LICENSE).
