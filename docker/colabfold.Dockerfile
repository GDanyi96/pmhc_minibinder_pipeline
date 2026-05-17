# ColabFold runtime — pinned CUDA-aware image.
# Pulled by bootstrap.sh.
FROM colabfold/colabfold:1.5.5-cuda12.2.2

WORKDIR /workspace
# AF2-multimer config (num_recycles=6 per MARES_IOANNIDIS_2025) is passed at runtime.
