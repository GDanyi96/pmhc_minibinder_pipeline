# pMHC-I Minibinder Pipeline — State

One-screen "where are we." Rulebook + metric/coding conventions + locked design
decisions: **`CLAUDE.md`**. Trap log: **`docs/known_traps.md`**. Portfolio (the
deliverable reviewers read): **`docs/`**. This file holds only current status
plus a few non-drifting reference tables that have no other home.

*Last updated: 2026-05-29 — cycle 03 complete + merged to `main`.*

---

## Now

- **Repo**: `github.com/GDanyi96/pmhc_minibinder_pipeline`. Cycle 03 merged to `main`.
- **Compute**: RunPod **A100-SXM4-80GB** pod `925f65a88d93`, **US-WA-1**, network volume at `/workspace`. Never push from pod.
- **Cycle 1** — controls validated (full target).
- **Cycle 2** — complete. Hero `design_2079_seq00` (99aa 4HB) is **peptide-BLIND** (`iface_pep` ∞; binder 28–40 Å from peptide) — a confident MHC-framework binder, not a reader.
- **Cycle 3** — complete + merged. 152 BAKER scaffolds → 72/150 geometry-pass → 288 AF2 folds → **91/288 (32%) engage peptide → 1 control-grade reader** `design_3010_seq00` (`iface_pep` 2.10; reads WT1 N5+Y8 via binder R55). Framework bias is **structural** (charge hypothesis tested + falsified). Best overall interface (`design_3084_seq02`, ipTM 0.89) and the framework champions are not peptide readers.
- **Cycle 4** — not started.

---

## Pod reference (does not drift)

| Env | Location | Manager | Use | Critical pins |
|---|---|---|---|---|
| Main pipeline | `/workspace/pipeline/.venv` | uv | ColabFold, ProteinMPNN, Snakemake, AF2 fan-in, metrics | JAX 0.4.34; torch 2.x+cu130; **nvidia-cudnn-cu12==9.1.0.70** (matched to driver 550) |
| RFdiffusion | `/workspace/miniconda3/envs/SE3nv` | conda + pip overlay | RFdiffusion inference only | Python 3.9; **torch 1.9.0+cu111**; dgl 0.9.1; e3nn 0.3.5. **Every pip install MUST use `--no-deps`** (Trap #36) |

- **Before any pod work**: `export TMPDIR=/workspace/.cache/tmp` (30 GB overlay fills otherwise, Trap #39). `uv run python` for pipeline; SE3nv python only for RFdiffusion — don't cross.
- `/workspace/RFdiffusion` @ `2d0c003df46b9db41d119321f15403dec3716cd9`; `Complex_base_ckpt.pt` sha256 `76e4e260…3250bc`.

---

## Chain layouts (empirically locked — Traps #29/#31/#32)

| Context | Layout |
|---|---|
| Cycle 02 full target, de novo (RFdiffusion out) | A=HC, B=β2m, C=peptide, **D=binder** |
| Cycle 03 sub-run A truncated, partial diffusion (RFdiffusion out) | **A=binder**, B=HLA[1:180], C=peptide (scaffold IDs preserved) |
| Truncated controls / Stage 2 FASTA | `HLA:peptide:binder` → AF2 out A=HLA, B=peptide, C=binder |
| `LAYOUT_CHAINS` | full: `binder=D, peptide=(C,), mhc=(A,B)` · truncated: `binder=C, peptide=(B,), mhc=(A,)` |

Partial-diffusion contig: `[N-N/0 B1-180/0 C1-9]` — bare `N-N` = redesigned binder; letter-prefixed = preserved motif (Trap #32).

---

## Metric convention (post cycle-03 audit)

Rank on **interface-8 Å iPAE** (`min(PAE_ij,PAE_ji)` over binder↔target pairs ≤ 8 Å); decompose into `iface_pep` / `iface_mhc`. **`iface_pep` is the specificity axis** (`∞` = peptide-blind). Bennett position-slice (`ipae_slice*`) kept for baseline only. Full- vs truncated-target values are **not** comparable. Full detail + control bands in `CLAUDE.md`.
