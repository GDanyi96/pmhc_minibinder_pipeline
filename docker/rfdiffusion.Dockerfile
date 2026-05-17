# RFdiffusion runtime image — thin recipe over the upstream Rosetta image.
# Pulled by bootstrap.sh; we extend only for project-specific I/O conventions.
FROM rosettacommons/rfdiffusion:latest

WORKDIR /workspace
# Project I/O conventions are mounted at runtime; nothing baked in here.
