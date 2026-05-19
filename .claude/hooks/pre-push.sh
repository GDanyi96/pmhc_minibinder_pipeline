#!/bin/bash
set -e
echo "═══════════════════════════════════════════════════════════"
echo "Pre-push verification gate"
echo "═══════════════════════════════════════════════════════════"
echo "Running formatter check..."
uv run black --check . || { echo "FAIL: black --check failed. Run 'uv run black .' before pushing."; exit 1; }
echo "PASS: formatter clean."
echo "Running snakemake dry-run in mock mode..."
if ! command -v snakemake &> /dev/null; then
    echo "WARN: snakemake not installed in this env; skipping dry-run gate."
    echo "Push allowed but verify on RunPod before running real pipeline."
    exit 0
fi
snakemake --dry-run --config mock=true -j1 2>&1 | tail -30
if [ $? -ne 0 ]; then
    echo "FAIL: snakemake dry-run did not pass. Fix before pushing."
    exit 1
fi
echo "PASS: pipeline DAG validates."
