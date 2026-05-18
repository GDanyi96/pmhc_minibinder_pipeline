# Stage 1 — RFdiffusion Backbone Generation

## Goal

Generate ~10 000 de novo minibinder backbones against the cleaned target,
biased by hotspots and constrained in length. For cycles ≥ 2, repurpose
previous winners via partial diffusion. Output is a directory of binder-only
backbone PDBs, ready for sequence design in stage 2.

## Inputs

- `data/targets/3hpj_clean.pdb` (from stage 0).
- `results/cycle_NN/stage0/target.yaml` — resolved hotspots, binder length range.
- `configs/rfdiffusion_default.yaml` — noise scaling, num_designs, contig template.
- (Cycle ≥ 2 only) `results/cycle_{NN-1}/winners/*.pdb` — previous-cycle survivors.

## Outputs

- `results/cycle_NN/stage1/backbones/design_{i:05d}.pdb` — N backbone PDBs
  (target + binder, but only binder coordinates are RFdiffusion output).
- `results/cycle_NN/stage1/manifest.json` — list of designs, hotspot used,
  diffusion seed, contig.

## Tools

- RFdiffusion installed natively on the pod (no DinD); upstream
  `rosettacommons/RFdiffusion` repo cloned by `bootstrap.sh`.
- Model weights: `/workspace/models/rfdiffusion/Complex_base_ckpt.pt`
  (downloaded by `bootstrap.sh`).
- Snakemake checkpoint to fan out per-design downstream steps.

## Anchored references

- `HADRUP_JENKINS_2025` — primary protocol: `noise_scale_ca=0`, hotspot-guided,
  contig schema `<target>/0 <binder_len>-<binder_len>`.
- `BAKER_LAB_2025` — partial diffusion for scaffold repurposing (`partial_T`
  noise re-noising of prior winners) across 11 pMHC-Is.

## Implementation tasks

- [ ] `run_rfdiffusion.py`: build the contig string from `target.yaml` +
      thresholds (binder length range).
- [ ] Compose the RFdiffusion CLI: `inference.num_designs`,
      `denoiser.noise_scale_ca=0`, `ppi.hotspot_res=[…]`,
      `contigmap.contigs=[…]`.
- [ ] Mount target PDB read-only into the container; mount
      `results/cycle_NN/stage1/` writable.
- [ ] Stream outputs into per-design PDB filenames; emit manifest.json.
- [ ] Cycle ≥ 2: switch to partial diffusion mode — feed
      `starting_pdbs_glob`, set `partial_T` from config.
- [ ] `--mock` flag: copy `tests/fixtures/rfdiffusion/sample.pdb` once to
      `design_00000.pdb`.

## Verification criteria

- Smoke test: `python workflow/scripts/run_rfdiffusion.py --config
  results/cycle_01/stage0/target.yaml --out
  results/cycle_01/stage1/ --mock` exits 0; produces one stub PDB.
- Real test (pod): GPU utilisation visible in `nvidia-smi` during run;
  10 000 backbones complete in <12 h on a single A100 80GB.

## Pitfalls

- RFdiffusion expects hotspots as `[A65, A150]` literal string in CLI; YAML
  list must be serialised carefully.
- Default RFdiffusion checkpoint is monomer; **must** use `Complex_base_ckpt`
  for pMHC binder design.
- `noise_scale_ca=0` is non-default. Confirm logs echo back the override.
- Output PDBs contain only binder coordinates by default; downstream
  ProteinMPNN expects the target re-attached. Stage 2a handles re-attachment.
- Snakemake checkpoint pattern required because N is only known after the
  rule runs (`num_designs` can be capped by --config override).
