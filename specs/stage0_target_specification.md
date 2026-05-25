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

## Sub-run A target: BAKER chain layout truncation

Cycle 03 sub-run A reuses BAKER's published pMHC-I binder scaffold library
(`BAKER_LAB_2025`). Each scaffold carries the target as a single **chain B =
HLA-A\*02:01 α1+α2 (residues 1-180) fused directly to the 9-mer peptide**
(189 continuous residues), with the binder on chain A — there is no β2m and no
α3 domain. Aligning such a scaffold onto our canonical `3hpj_clean.pdb`
(chain B = β2m) superposes the HLA groove onto β2m, two unrelated folds, and
silently yields garbage starting geometry (see `docs/known_traps.md` Trap #31).

### Output

`prep_baker_target.py` derives `data/targets/3hpj_baker_truncated.pdb` from
`3hpj_clean.pdb`:

- keep chain A residues 1-180 (the α1+α2 peptide-binding groove); drop 181-275
  (α3),
- rename that chain A → **B** (BAKER convention),
- drop the original chain B (β2m) entirely,
- keep chain C (peptide RMFPNAPYL, residues 1-9) verbatim,
- ATOM records only (HETATM/water/H/alt-loc-B stripped); original numbering
  preserved (no renumber). Validated CA counts: **180 on B, 9 on C**.

Sub-run B (`design_2079`) keeps the full `3hpj_clean.pdb`: it was AF2-validated
against the full 4-chain target and must not change.

### Geometric proof (binding groove preserved)

The pMHC-I peptide-binding cleft is formed entirely by the HLA α1 (residues
~1-90) and α2 (~91-180) domains, which present the floor β-sheet and the two
α-helices that wall the peptide. β2m and α3 are membrane-proximal Ig domains
that do **not** contact the bound peptide and do not line the groove. Removing
them therefore leaves every binder-relevant atom — the entire α1/α2 groove plus
the peptide — in place, so a binder designed against the truncated fragment
faces an identical interface to one designed against the full target. The
peptide-facing hotspots used by sub-run A (`C1,C4,C6,C8` on the peptide;
`B65,B66,B150,B155` on the α-helix walls after A→B rename) all fall within the
retained 1-180 span.

### Chain-layout difference between the two paths

| Path                 | Reference / fold target              | AF2 output chains                         | iPAE chains (LAYOUT_CHAINS)                          |
|----------------------|--------------------------------------|-------------------------------------------|------------------------------------------------------|
| **full** (default)   | `3hpj_clean.pdb` (4-chain)           | A=HC, B=β2m, C=peptide, D=binder          | binder=D, peptide=(C,), mhc=(A,B), combined=(A,B,C)  |
| **truncated** (sub-run A) | `3hpj_baker_truncated.pdb` (B=HLA, C=peptide) | A=HLA[1:180], B=peptide, C=binder | binder=C, peptide=(B,), mhc=(A,), combined=(A,B)     |

Truncated controls (`run_controls.py --target truncated`) fold the binder
against the β2m-free groove (`HLA : peptide : binder` FASTA), matching sub-run
A's starting geometry, and write to
`results/cycle_NN/controls_truncated_baseline/metrics_truncated_baseline.json`.
`compute_metrics` selects the binder/peptide/MHC chains via
`LAYOUT_CHAINS[target_layout]`, so decomposed iPAE
(`ppi_pae_int_peptide` / `ppi_pae_int_mhc`) is computed correctly for both
layouts with no behavior change on the full path.

### Anchored references

- `BAKER_LAB_2025` — scaffold library chain layout (binder + fused HLA-peptide
  target chain); upstream `align_chainB.py` methodology.
