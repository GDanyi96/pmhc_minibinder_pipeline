#!/usr/bin/env python3
"""design_3010_seq00 — per-peptide-residue contact + AF2 confidence.

RMFPNAPYL on HLA-A*02:01: P2-Met and P9-Leu are buried anchors; the central
bulge (P4-P8) is solvent/TCR-facing and is what a peptide-SPECIFIC binder must
read. Reports min heavy-atom distance binder->each peptide residue, the closest
binder residue, and the mean PAE for that peptide residue vs the binder.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser

ROOT = Path("/workspace/pipeline")
PRED = ROOT / "results/cycle_03/stage2/subrun_a/af2_designs/predictions"
SEQ = "design_3010_seq00"
ANNOT = {1:("R","N-term"),2:("M","ANCHOR buried"),3:("F","semi-exposed"),
         4:("P","EXPOSED"),5:("N","EXPOSED-specificity"),6:("A","EXPOSED"),
         7:("P","EXPOSED"),8:("Y","EXPOSED-specificity"),9:("L","ANCHOR buried")}

sj  = sorted(PRED.glob(f"**/{SEQ}_scores_rank_001*.json"))[0]
pdb = sorted(PRED.glob(f"**/{SEQ}_unrelaxed_rank_001*.pdb"))[0]
print(f"PDB: {pdb.relative_to(ROOT)}")
s = PDBParser(QUIET=True).get_structure("d", str(pdb))[0]

# confirm chain layout + peptide identity
from Bio.SeqUtils import seq1
def chain_seq(c): return "".join(seq1(r.get_resname()) for r in s[c] if r.id[0].strip()=="")
pep_seq = chain_seq("B")
print(f"chain A len={sum(1 for r in s['A'] if r.id[0].strip()=='')} (HLA)")
print(f"chain B seq={pep_seq} (peptide; expect RMFPNAPYL)")
print(f"chain C len={sum(1 for r in s['C'] if r.id[0].strip()=='')} (binder)")
assert pep_seq == "RMFPNAPYL", f"peptide chain mismatch: got {pep_seq}"

pep = [(r.id[1], [a for a in r if a.element!='H']) for r in s["B"] if r.id[0].strip()==""]
bind_atoms = [(r.id[1], a) for r in s["C"] if r.id[0].strip()=="" for a in r if a.element!='H']
bcoords = np.array([a.get_coord() for _,a in bind_atoms]); bres = [rn for rn,_ in bind_atoms]

pae = np.array(json.loads(sj.read_text())["pae"])
nA = sum(1 for r in s['A'] if r.id[0].strip()=="")
nB = len(pep); nC = len(bres) and sum(1 for r in s['C'] if r.id[0].strip()=="")
cs, ce = nA+nB, nA+nB+nC

print(f"\n{'P':>2} {'aa':>2} {'role':<20} {'min_dist':>9} {'closest_binder_res':>18} {'PAE_to_binder':>13} {'contact':>8}")
print("-"*82)
contacted_exposed, contacted_anchor = [], []
for i,(rn,atoms) in enumerate(pep, start=1):
    coords = np.array([a.get_coord() for a in atoms])
    d = np.linalg.norm(coords[:,None,:]-bcoords[None,:,:],axis=-1)
    dmin = d.min(); cb = bres[d.min(axis=0).argmin()]
    pidx = nA + (i-1)
    paeb = (pae[pidx, cs:ce].mean()+pae[cs:ce, pidx].mean())/2
    aa,role = ANNOT[i]
    tag = "TIGHT" if dmin<=5 else ("contact" if dmin<=8 else "")
    print(f"{i:>2} {aa:>2} {role:<20} {dmin:>9.2f} {cb:>18} {paeb:>13.2f} {tag:>8}")
    if dmin<=8 and "EXPOSED" in role: contacted_exposed.append(f"{aa}{i}")
    if dmin<=8 and "ANCHOR" in role: contacted_anchor.append(f"{aa}{i}")

print("\nVERDICT")
print(f"  exposed/specificity residues contacted (<=8A): {contacted_exposed or 'NONE'}")
print(f"  anchor residues contacted (<=8A):              {contacted_anchor or 'NONE'}")
print(f"  -> genuinely WT1-specific if it engages the central bulge (P4-P8), esp. N5/Y8")

# binder residues forming the interface (for the eventual design/mutagenesis story)
near = set()
for _,atoms in pep:
    coords = np.array([a.get_coord() for a in atoms])
    d = np.linalg.norm(coords[:,None,:]-bcoords[None,:,:],axis=-1)
    for j in np.where(d.min(axis=0)<=5)[0]: near.add(bres[j])
print(f"  binder residues within 5A of peptide: {sorted(near)}")
