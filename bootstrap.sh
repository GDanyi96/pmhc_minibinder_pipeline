#!/usr/bin/env bash
# One-command RunPod setup for the pMHC-I minibinder pipeline.
# Usage (on the pod):  cd /workspace/pipeline && bash bootstrap.sh
#   --with-rfdiffusion-weights  also fetch Stage 1 weights (defer until Stage 1)

set -euo pipefail

WITH_RFDIFF_WEIGHTS=0
for arg in "$@"; do
  case "$arg" in
    --with-rfdiffusion-weights) WITH_RFDIFF_WEIGHTS=1 ;;
    -h|--help)
      sed -n '2,4p' "$0"
      exit 0
      ;;
    *) echo "bootstrap.sh: unknown arg '$arg'" >&2; exit 2 ;;
  esac
done

# RunPod base sets LD_LIBRARY_PATH=/usr/local/cuda/lib64 (CUDA 13). JAX
# cu12 wheels bundle their own CUDA 12 libs; the system override loads
# mismatched libs into the JAX process → SIGSEGV in cuDNN at first
# forward pass. Unset for this shell and all child processes (uv run,
# colabfold_batch subprocess). Empirically reproduced on cycle 1 pod run.
unset LD_LIBRARY_PATH

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
echo "[3/7] Installing Python deps via uv (all extras: dev + colabfold + proteinmpnn)..."
uv sync --all-extras

# ─── 4. Native tool installs (no Docker — the pod IS the container) ─
echo
echo "[4/7] Setting up native tool installs..."
echo "  The pod runs all heavy tools natively (no Docker-in-Docker)."
echo "  ColabFold (colabfold_batch CLI) was installed via the 'colabfold' extra in step 3."

# ProteinMPNN ships as a git repo, not a pip wheel. Clone (or update) under
# /workspace so run_proteinmpnn.py's default PROTEINMPNN_DIR resolves.
# The vanilla_model_weights/ subdir ships with the clone (no separate
# download URL); we sha256-verify the bundled weight file in-place.
PROTEINMPNN_DIR="${PROTEINMPNN_DIR:-/workspace/ProteinMPNN}"
if [ -d "$PROTEINMPNN_DIR/.git" ]; then
  echo "  ProteinMPNN already cloned at $PROTEINMPNN_DIR — pulling latest..."
  git -C "$PROTEINMPNN_DIR" pull --ff-only
else
  echo "  Cloning ProteinMPNN to $PROTEINMPNN_DIR..."
  git clone https://github.com/dauparas/ProteinMPNN "$PROTEINMPNN_DIR"
fi

# Pinned sha256 for ProteinMPNN's bundled v_48_020.pt weight file.
# Self-bootstrapping: leave empty in the initial commit; on first run the
# script prints the computed value and aborts asking the operator to pin
# it. Same flow applies when upstream rotates the weight.
PROTEINMPNN_CKPT_SHA256="${PROTEINMPNN_CKPT_SHA256:-}"
verify_proteinmpnn_weight() {
  local ckpt="$PROTEINMPNN_DIR/vanilla_model_weights/v_48_020.pt"
  if [ ! -f "$ckpt" ]; then
    echo "FAIL: $ckpt not found — did the clone include vanilla_model_weights/?"
    exit 1
  fi
  if [ -z "$PROTEINMPNN_CKPT_SHA256" ]; then
    local computed
    computed=$(sha256sum "$ckpt" | cut -d' ' -f1)
    echo "FAIL: PROTEINMPNN_CKPT_SHA256 not pinned in bootstrap.sh."
    echo "Computed value for this checkout: $computed"
    echo "Edit bootstrap.sh and set PROTEINMPNN_CKPT_SHA256=\"$computed\", then re-run."
    exit 1
  fi
  if ! echo "${PROTEINMPNN_CKPT_SHA256}  ${ckpt}" | sha256sum -c - >/dev/null; then
    echo "FAIL: v_48_020.pt sha256 mismatch — upstream rotated the weight?"
    echo "  expected: $PROTEINMPNN_CKPT_SHA256"
    echo "  actual:   $(sha256sum "$ckpt" | cut -d' ' -f1)"
    exit 1
  fi
  echo "  v_48_020.pt: sha256 verified"
}
verify_proteinmpnn_weight

# Sanity: colabfold_batch must be on PATH after the colabfold extra install.
if uv run colabfold_batch --help >/dev/null 2>&1; then
  echo "  colabfold_batch: OK"
else
  echo "  WARN: colabfold_batch --help failed; check the colabfold extra install."
fi

# ─── 5. Native tool installs (Stage 1) + model weights (opt-in) ─
echo
echo "[5/7] RFdiffusion clone + weights..."

# RFdiffusion ships as a git repo, not a pip wheel. Clone (or update) under
# /workspace so run_rfdiffusion.py's hardcoded /workspace/RFdiffusion path
# resolves. The clone is unconditional (cheap; ~50 MB), the weights are
# gated by --with-rfdiffusion-weights because the Complex_base_ckpt.pt
# weight file is ~1.5 GB.
clone_rfdiffusion() {
  local dir="${RFDIFFUSION_DIR:-/workspace/RFdiffusion}"
  if [ -d "$dir/.git" ]; then
    echo "  RFdiffusion already cloned at $dir — pulling latest..."
    git -C "$dir" pull --ff-only
  else
    echo "  Cloning RFdiffusion to $dir..."
    git clone https://github.com/RosettaCommons/RFdiffusion "$dir"
  fi
}
clone_rfdiffusion

# Pinned sha256 for Complex_base_ckpt.pt. Self-bootstrapping: leave empty
# in the initial commit; on first download the script prints the computed
# value and aborts asking the operator to pin it. Same flow applies in
# cycle 3 if the upstream checkpoint rotates.
RFDIFF_CKPT_SHA256="${RFDIFF_CKPT_SHA256:-}"
# URL scheme note: RFdiffusion README publishes http://files.ipd.uw.edu/...;
# verify whether https://files.ipd.uw.edu/... resolves before final pin.
# sha256 verification closes any in-transit integrity gap regardless.
RFDIFF_CKPT_URL="${RFDIFF_CKPT_URL:-http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Complex_base_ckpt.pt}"

download_rfdiffusion_weights() {
  local models_dir="${RFDIFFUSION_MODELS_DIR:-/workspace/models/RFdiffusion}"
  mkdir -p "$models_dir"
  local ckpt="$models_dir/Complex_base_ckpt.pt"

  if [ -f "$ckpt" ]; then
    if [ -n "$RFDIFF_CKPT_SHA256" ]; then
      if echo "${RFDIFF_CKPT_SHA256}  ${ckpt}" | sha256sum -c - >/dev/null 2>&1; then
        echo "  Complex_base_ckpt.pt: already present, sha256 OK"
        return 0
      fi
      echo "  Complex_base_ckpt.pt: already present but sha256 mismatch — re-downloading"
      rm -f "$ckpt"
    else
      echo "  Complex_base_ckpt.pt: already present, no pinned sha256 to verify"
      return 0
    fi
  fi

  echo "  Downloading Complex_base_ckpt.pt from $RFDIFF_CKPT_URL..."
  wget -q "$RFDIFF_CKPT_URL" -O "$ckpt"

  if [ -z "$RFDIFF_CKPT_SHA256" ]; then
    local computed
    computed=$(sha256sum "$ckpt" | cut -d' ' -f1)
    echo "FAIL: RFDIFF_CKPT_SHA256 not pinned in bootstrap.sh."
    echo "Computed value for this download: $computed"
    echo "Edit bootstrap.sh and set RFDIFF_CKPT_SHA256=\"$computed\", then re-run."
    exit 1
  fi

  if ! echo "${RFDIFF_CKPT_SHA256}  ${ckpt}" | sha256sum -c -; then
    echo "FAIL: Complex_base_ckpt.pt sha256 mismatch — abort"
    rm -f "$ckpt"
    exit 1
  fi
  echo "  Complex_base_ckpt.pt: downloaded, sha256 verified"
}

if [ "$WITH_RFDIFF_WEIGHTS" -eq 1 ]; then
  download_rfdiffusion_weights
  mkdir -p /workspace/models/esm2
  echo "  TODO: huggingface-cli download facebook/esm2_t33_650M_UR50D -> /workspace/models/esm2/"
else
  echo "  Stage 1 weights + Stage 4 weights skipped (use --with-rfdiffusion-weights to fetch)."
fi

# ─── 6. Target PDBs (download + clean) ─────────────────────────
echo
echo "[6/7] Target PDBs..."
mkdir -p data/targets
for pdb in 3HPJ 2BNR; do
  raw="data/targets/${pdb}.pdb"
  if [ -f "$raw" ]; then
    echo "  $pdb: already present, skipping download"
  else
    wget -q "https://files.rcsb.org/download/${pdb}.pdb" -O "$raw"
    echo "  $pdb: downloaded"
  fi
done

clean_primary="data/targets/3hpj_clean.pdb"
clean_control="data/targets/2bnr_clean.pdb"
target_yaml="data/targets/target.yaml"
if [ -f "$clean_primary" ] && [ -f "$clean_control" ] && [ -f "$target_yaml" ]; then
  echo "  cleaned PDBs + target.yaml present, skipping prep_target"
else
  echo "  Running prep_target.py..."
  uv run python workflow/scripts/prep_target.py \
    --primary-config configs/target_wt1_a0201.yaml \
    --control-config configs/target_2bnr_a0201.yaml \
    --primary-pdb-in data/targets/3HPJ.pdb \
    --control-pdb-in data/targets/2BNR.pdb \
    --primary-pdb-out "$clean_primary" \
    --control-pdb-out "$clean_control" \
    --out "$target_yaml"
fi

# ─── 7. Smoke test ─────────────────────────────────────────────
echo
echo "[7/7] Running snakemake dry-run in mock mode..."
uv run snakemake --dry-run --config mock=true -j1 2>&1 | tail -5

echo
echo "═══════════════════════════════════════════════════════════"
echo "READY"
echo "═══════════════════════════════════════════════════════════"
echo "Next:"
echo "  uv run python scripts/run_controls.py   # Stage 2 controls (~10 min on A100)"
echo "  uv run snakemake --cores all            # full cycle"
