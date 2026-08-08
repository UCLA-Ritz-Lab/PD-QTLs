#!/usr/bin/env bash
# Run the 1000G-anchored ancestry PC projection end to end.
#
# Safe to re-run: Snakemake picks up where it left off, and the per-chromosome
# 1000G VCFs are temp() so they are removed once merged rather than being
# re-downloaded. Logs land in ../../../results/ANCESTRY_REF/logs/.
set -euo pipefail

cd "$(dirname "$0")"
# shellcheck disable=SC1091
source /home/garyc/miniconda3/etc/profile.d/conda.sh
conda activate pd-qtl

exec snakemake \
  --sdm conda \
  --cores "${CORES:-4}" \
  --rerun-incomplete \
  --printshellcmds \
  "$@"
