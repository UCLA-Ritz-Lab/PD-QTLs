#!/usr/bin/env bash
# Entrypoint for each Google Batch task. Google Batch sets BATCH_TASK_INDEX
# (0-based) automatically for every task in the job. We use it to pick one
# line out of manifest.txt -- one Sentrix_ID per line -- so each task
# processes exactly one batch directory with process_batch.R.
set -euo pipefail

: "${BUCKET:?BUCKET env var required}"
: "${RAW_PREFIX:?RAW_PREFIX env var required}"
: "${NORM_PREFIX:?NORM_PREFIX env var required}"
: "${BETAS_PREFIX:?BETAS_PREFIX env var required}"
: "${BATCH_TASK_INDEX:?Not running under Google Batch -- BATCH_TASK_INDEX unset}"
ARRAY_TYPE="${ARRAY_TYPE:-EPIC}"

MANIFEST_LOCAL=/tmp/manifest.txt
gcloud storage cp "gs://${BUCKET}/${RAW_PREFIX}/manifest.txt" "$MANIFEST_LOCAL"

# BATCH_TASK_INDEX is 0-based; sed is 1-based.
LINE_NUM=$((BATCH_TASK_INDEX + 1))
SENTRIX_ID=$(sed -n "${LINE_NUM}p" "$MANIFEST_LOCAL")

if [[ -z "$SENTRIX_ID" ]]; then
  echo "No Sentrix_ID found for task index ${BATCH_TASK_INDEX} (line ${LINE_NUM}) -- exiting."
  exit 1
fi

echo "=== Task ${BATCH_TASK_INDEX}: processing Sentrix batch ${SENTRIX_ID} ==="

# process_batch.R does file.path(opt$idat_dir, opt$sentrix_id) internally,
# so idat_dir must be the PARENT of the Sentrix_ID folder, not the folder
# itself. out_dir is a shared results root -- the script namespaces every
# output file by sentrix_id, so multiple tasks writing into the same local
# root (and later the same gs:// prefix) don't collide.
LOCAL_IDAT_PARENT="/tmp/idats"
LOCAL_OUT="/tmp/normalized"
mkdir -p "${LOCAL_IDAT_PARENT}/${SENTRIX_ID}" "$LOCAL_OUT"

gcloud storage cp -r "gs://${BUCKET}/${RAW_PREFIX}/${SENTRIX_ID}/*" "${LOCAL_IDAT_PARENT}/${SENTRIX_ID}/"

# process_batch.R opens the --log file directly without creating its parent
# dir, so create it ahead of time.
mkdir -p "${LOCAL_OUT}/logs"

Rscript /pipeline/scripts/process_batch.R \
  --sentrix_id "$SENTRIX_ID" \
  --idat_dir "$LOCAL_IDAT_PARENT" \
  --out_dir "$LOCAL_OUT" \
  --array_type "$ARRAY_TYPE" \
  --log "${LOCAL_OUT}/logs/${SENTRIX_ID}.log"

gcloud storage cp -r "${LOCAL_OUT}/*" "gs://${BUCKET}/${NORM_PREFIX}/"

# Convert this batch's betas RDS -> CSV.gz inline, right here in the
# fan-out step, rather than as a separate sequential Snakemake rule on the
# downstream VM. This step is per-batch/embarrassingly parallel just like
# process_batch.R itself, so there's no reason to pay for a second Batch
# job (or a slow sequential loop on the downstream VM) to do it.
RDS_PATH="${LOCAL_OUT}/betas/${SENTRIX_ID}.betas.rds"
CSV_PATH="${LOCAL_OUT}/betas_csv/${SENTRIX_ID}.betas.csv.gz"
mkdir -p "${LOCAL_OUT}/betas_csv"

Rscript /pipeline/scripts/convert_rds_to_csv.R \
  --rds "$RDS_PATH" \
  --csv "$CSV_PATH"

gcloud storage cp "$CSV_PATH" "gs://${BUCKET}/${BETAS_PREFIX}/${SENTRIX_ID}.betas.csv.gz"

echo "=== Task ${BATCH_TASK_INDEX}: done, uploaded ${SENTRIX_ID} results (RDS + converted CSV) ==="
