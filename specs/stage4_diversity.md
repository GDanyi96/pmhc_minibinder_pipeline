# Stage 4 — Diversity Selection (ESM-2 + Farthest-Point Sampling)

## Goal

Embed every specificity-passing design with ESM-2 (650M), mean-pool over
binder residues, and select the top-N most diverse designs via
farthest-point sampling (FPS). Compare coverage to a k-means baseline.
FPS is chosen over k-means because *coverage* of sequence space matters
more than density when designs are cheap to generate but experimental
slots are scarce.

## Inputs

- `results/cycle_NN/stage3/specificity_survivors.json` — survivors of
  cross-pan filtering.
- `results/cycle_NN/stage2/proteinmpnn/design_*.fasta` — binder sequences.
- `configs/thresholds.yaml` (`diversity.fps_n_select`) — N to keep (96).

## Outputs

- `results/cycle_NN/stage4/embeddings.npz` — `{ids: [...], embeddings:
  (N, 1280)}`.
- `results/cycle_NN/stage4/fps_selection.json` — IDs of selected designs.
- `results/cycle_NN/stage4/diversity_report.json` — coverage radius,
  mean pairwise distance, comparison vs k-means baseline.

## Tools

- `fair-esm` (`esm.pretrained.esm2_t33_650M_UR50D`).
- `torch` (CUDA).
- `scikit-learn` for k-means baseline.
- `numpy` for FPS (implementation < 50 lines).

## Anchored references

- `HOUSEHOLDER_GARCIA_2025` — uses ESM-derived embeddings to filter binder
  candidates for proteome-wide off-target safety; we use the same model
  but for forward diversity selection.

## Implementation tasks

- [ ] `embed_designs.py`: load ESM-2 650M, batch sequences (16 at a time
      on A100), extract layer-33 representations, mean-pool over binder
      residue positions.
- [ ] Save embeddings.npz (float16 to save space).
- [ ] FPS: pick farthest design from random seed, iteratively add the
      design with max-min Euclidean distance to selected set, until N.
- [ ] K-means baseline: cluster into N groups, take cluster centroids;
      compute coverage radius vs FPS.
- [ ] Emit diversity_report.json with both metrics.
- [ ] `--mock`: cp `tests/fixtures/embeddings/sample.npz` to output.

## Verification criteria

- Smoke test: `embed_designs.py --mock` exits 0 in <1s.
- Real test (pod): 1 000 designs embedded in <5 min on A100.
- FPS coverage radius < k-means coverage radius (sanity check).

## Pitfalls

- ESM-2 650M takes ~3 GB VRAM; safe alongside RFdiffusion (which is idle
  during this stage). Don't run concurrently with stage 1.
- Mean-pool over **binder** residues only — sequences passed in are
  already binder-only, but assert length matches binder spec.
- FPS is sensitive to the random seed for the initial point. Use a
  deterministic seed (e.g., 42) and log it.
- Float16 embeddings: keep float32 inside FPS to avoid distance underflow;
  cast to float16 only on disk.
- K-means in 1280-D space is high-variance; run 10 random restarts and
  take the best inertia for a fair baseline.
