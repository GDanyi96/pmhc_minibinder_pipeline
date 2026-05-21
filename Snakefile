# Top-level pipeline DAG.
# Per-stage rule files live under workflow/rules/.
# Real-mode scripts live under workflow/scripts/.
# Mock-mode populates outputs from tests/fixtures/ — every rule supports
# `--config mock=true` for a no-GPU end-to-end smoke test.

from pathlib import Path

# Default config: real-mode, cycle 01.
config.setdefault("mock", False)
config.setdefault("cycle", "01")

CYCLE: str = str(config["cycle"]).zfill(2)
MOCK: bool = bool(config["mock"])
RESULTS = Path(f"results/cycle_{CYCLE}")
FIXTURES = Path("tests/fixtures")

# Stage 1 backbone-generation mode. Cycle 03 replaces de-novo RFdiffusion
# with partial diffusion from two seeded sub-runs (BAKER scaffold library +
# our cycle-02 hero), merged into the canonical Stage 1 output. Defaults from
# the cycle number; overridable via `--config stage1_mode=denovo|partial`.
config.setdefault("stage1_mode", "partial" if int(CYCLE) >= 3 else "denovo")
STAGE1_PARTIAL: bool = config["stage1_mode"] == "partial"

# Cycle 03 designs use the Ala-biased ProteinMPNN config (Trap #30);
# earlier cycles keep the unbiased default.
MPNN_DESIGNS_CONFIG: str = (
    "configs/proteinmpnn_cycle03.yaml"
    if int(CYCLE) >= 3
    else "configs/proteinmpnn_default.yaml"
)

include: "workflow/rules/00_target_prep.smk"
include: "workflow/rules/01_rfdiffusion.smk"
include: "workflow/rules/01_stage1_subrun_a.smk"
include: "workflow/rules/01_stage1_subrun_b.smk"
include: "workflow/rules/01_stage1_merge.smk"
include: "workflow/rules/02_proteinmpnn.smk"
include: "workflow/rules/03_colabfold.smk"
include: "workflow/rules/04_metrics.smk"
include: "workflow/rules/02b_proteinmpnn_designs.smk"
include: "workflow/rules/03b_af2_designs.smk"
include: "workflow/rules/05_crosspan.smk"
include: "workflow/rules/06_embedding.smk"
include: "workflow/rules/07_active_learning.smk"
include: "workflow/rules/08_reporting.smk"

# `stage1_merge` and `rfdiffusion` both produce stage1/rfdiffusion/
# {stage1_summary.json,designs.jsonl}. Resolve the ambiguity by mode: partial
# diffusion (cycle >= 3) uses the merge of sub-runs A+B; otherwise the legacy
# de-novo rule wins and the sub-run rules stay inert (never requested).
if STAGE1_PARTIAL:
    ruleorder: stage1_merge > rfdiffusion
else:
    ruleorder: rfdiffusion > stage1_merge


rule all:
    input:
        f"reports/cycle_{CYCLE}.md",
        str(RESULTS / "stage2" / ".done"),
