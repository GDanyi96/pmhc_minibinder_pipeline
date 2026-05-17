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

include: "workflow/rules/00_target_prep.smk"
include: "workflow/rules/01_rfdiffusion.smk"
include: "workflow/rules/02_proteinmpnn.smk"
include: "workflow/rules/03_colabfold.smk"
include: "workflow/rules/04_metrics.smk"
include: "workflow/rules/05_crosspan.smk"
include: "workflow/rules/06_embedding.smk"
include: "workflow/rules/07_active_learning.smk"
include: "workflow/rules/08_reporting.smk"


rule all:
    input:
        f"reports/cycle_{CYCLE}.md",
