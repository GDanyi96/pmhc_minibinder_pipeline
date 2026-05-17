# Stage 0 — Target Specification

## Goal

Resolve the primary target (WT1 RMFPNAPYL / HLA-A\*02:01, PDB 3HPJ) and the
positive-control target (Jenkins NY1-B04 vs SLLMWITQC, PDB 2BNR) into clean,
chain-renamed, hetatm-stripped PDB files with a confirmed hotspot list ready
for RFdiffusion. This stage is the input contract for stage 1.

## Inputs

- `configs/target_wt1_a0201.yaml` — target metadata (chains, hotspot
  candidates, binder length range).
- `data/targets/3HPJ.pdb` — raw RCSB download (fetched by `bootstrap.sh`).
- `data/targets/2BNR.pdb` — positive-control raw download.

## Outputs

- `data/targets/3hpj_clean.pdb` — chains renamed to {A: HC, B: β2m, C:
  peptide}, hetatms/waters stripped, single biological assembly.
- `data/targets/2bnr_clean.pdb` — same conventions.
- `results/cycle_01/stage0/target.yaml` — resolved hotspot residues
  (numbered per cleaned PDB), peptide sequence, binder chain ID.

## Tools

- `biopython` (`Bio.PDB`) for parsing + chain manipulation.
- `pydantic` for YAML schema validation.
- `pymol` (optional) to validate hotspot solvent accessibility — not required for v0.

## Anchored references

- `HADRUP_JENKINS_2025` — chain ordering convention (HC=A, β2m=B, peptide=C,
  binder=D); hotspot definition on peptide-facing α1/α2 helices.
- `BAKER_LAB_2025` — confirms binder chain is appended last with a >12 Å gap
  in residue numbering to avoid AF2 chain-merge artefacts.

## Implementation tasks

- [ ] `prep_target.py`: parse raw PDB, retain biological assembly 1.
- [ ] Rename chains to A=HC, B=β2m, C=peptide; verify lengths
      (HC ≈ 275, β2m ≈ 99, peptide = 9).
- [ ] Strip hetatms, waters, alt-loc B, hydrogens.
- [ ] Renumber peptide residues 1–N (preserve original in REMARK).
- [ ] Validate hotspot list from `configs/target_*.yaml` resolves to real
      residues in cleaned PDB; reject otherwise.
- [ ] Emit `results/cycle_NN/stage0/target.yaml` with resolved hotspots,
      sequences, paths.
- [ ] `--mock` flag: copy `tests/fixtures/target_3hpj_clean.pdb` to output.

## Verification criteria

- Smoke test: `python workflow/scripts/prep_target.py --config
  configs/target_wt1_a0201.yaml --out results/cycle_01/stage0/target.yaml
  --mock` exits 0 in <1s and writes a non-empty target.yaml stub.
- Real test (on pod): cleaned 3HPJ has exactly 3 chains, peptide chain has
  9 residues, all hotspot residues in `target.yaml` exist in cleaned PDB.

## Pitfalls

- 3HPJ has multiple biological assemblies — must pick assembly 1.
- `Bio.PDB` will silently rename chains on save if IDs collide. Verify
  chain IDs match expectation after save.
- AF2-multimer is sensitive to chain ordering; downstream stages assume the
  D-chain (binder) is appended last. Do not reorder.
- Hotspot residue numbers in published RFdiffusion examples reference RCSB
  numbering, not cleaned/renumbered. Track both.
