#!/usr/bin/env bash
# One-command RunPod setup for the pMHC-I minibinder pipeline.
# Usage (on the pod):  cd /workspace/pipeline && bash bootstrap.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

echo "═══════════════════════════════════════════════════════════"
echo "pMHC-I minibinder pipeline — RunPod bootstrap"
echo "═══════════════════════════════════════════════════════════"

# ─── 1. GPU sanity check ───────────────────────────────────────
echo
echo "[1/7] Verifying GPU..."
if ! command -v nvidia-smi &>/dev/null; then
  echo "FAIL: nvidia-smi not found — is this actually a GPU pod?"
  exit 1
fi
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ─── 2. Install uv ─────────────────────────────────────────────
echo
echo "[2/7] Installing uv if missing..."
if ! command -v uv &>/dev/null; then
  curl -LsSf https://astral.sh/uv/install.sh | sh
  export PATH="$HOME/.local/bin:$PATH"
fi
uv --version

# ─── 3. Python deps ────────────────────────────────────────────
echo
echo "[3/7] Installing Python deps via uv..."
uv sync --extra dev

# ─── 4. Docker images (real downloads — uncomment when ready) ──
echo
echo "[4/7] Pulling Docker images..."
echo "TODO: docker pull rosettacommons/rfdiffusion:latest"
echo "TODO: docker pull colabfold/colabfold:1.5.5-cuda12.2.2"

# ─── 5. Model weights ──────────────────────────────────────────
echo
echo "[5/7] Downloading model weights to /workspace/models/..."
mkdir -p /workspace/models/rfdiffusion /workspace/models/esm2
echo "TODO: wget RFdiffusion Complex_base_ckpt.pt → /workspace/models/rfdiffusion/"
echo "TODO: huggingface-cli download facebook/esm2_t33_650M_UR50D → /workspace/models/esm2/"

# ─── 6. Target PDBs ────────────────────────────────────────────
echo
echo "[6/7] Downloading target PDBs..."
mkdir -p data/targets
if [ ! -f data/targets/3HPJ.pdb ]; then
  echo "TODO: wget https://files.rcsb.org/download/3HPJ.pdb -O data/targets/3HPJ.pdb"
fi
if [ ! -f data/targets/2BNR.pdb ]; then
  echo "TODO: wget https://files.rcsb.org/download/2BNR.pdb -O data/targets/2BNR.pdb"
fi

# ─── 7. Smoke test ─────────────────────────────────────────────
echo
echo "[7/7] Running snakemake dry-run in mock mode..."
uv run snakemake --dry-run --config mock=true -j1 2>&1 | tail -5

echo
echo "═══════════════════════════════════════════════════════════"
echo "READY"
echo "═══════════════════════════════════════════════════════════"
echo "Next: replace TODO lines above with real downloads, then:"
echo "  uv run snakemake --cores all"
