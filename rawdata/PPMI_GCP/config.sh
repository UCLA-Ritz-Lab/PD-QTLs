#!/usr/bin/env bash
# Shared config — edit these values, then `source config.sh` (all other
# scripts do this automatically).

# --- Project / region -------------------------------------------------
export PROJECT_ID="ppmi-batch-pipeline"        # suggest a fresh project, separate from pioneer-metals-501918
export REGION="us-west1"                       # keep consistent with Pioneer Metals for familiarity
export ZONE="us-west1-a"

# --- Storage -------------------------------------------------------
export BUCKET="ppmi-idat-pipeline-${PROJECT_ID}"   # must be globally unique
export RAW_PREFIX="raw/idats"                       # gs://$BUCKET/raw/idats/{SENTRIX_ID}/...
export NORM_PREFIX="results/normalized"             # gs://$BUCKET/results/normalized/{SENTRIX_ID}/...
export BETAS_PREFIX="results/betas"                 # gs://$BUCKET/results/betas/{sentrix_id}.betas.csv.gz -- matches local OUT/betas/ convention, produced inline in the fan-out step now instead of a separate convert_rds_to_csv Snakemake rule
export FINAL_PREFIX="results/final"                 # gs://$BUCKET/results/final/...

# --- Artifact Registry / container image ----------------------------
export AR_REPO="ppmi-pipeline"
export IMAGE_NAME="sesame-normalize"
export IMAGE_TAG="latest"
export IMAGE_URI="${REGION}-docker.pkg.dev/${PROJECT_ID}/${AR_REPO}/${IMAGE_NAME}:${IMAGE_TAG}"

# --- Service accounts -------------------------------------------------
export BATCH_SA="batch-runner@${PROJECT_ID}.iam.gserviceaccount.com"
export VM_SA="downstream-vm@${PROJECT_ID}.iam.gserviceaccount.com"

# --- Batch job sizing --------------------------------------------------
export NORMALIZE_CPU_MILLI=2000     # 2 vCPU per task — SeSAMe per-batch normalization is single-threaded-ish
export NORMALIZE_MEM_MIB=8192       # 8 GB per task
export NORMALIZE_MAX_PARALLEL=8     # cap concurrent Spot tasks (tune to your batch count / quota)
export ARRAY_TYPE="EPIC"            # PPMI methylation array is EPIC v1, not v2 — passed through to process_batch.R --array_type

# --- Downstream VM sizing -----------------------------------------------
export DOWNSTREAM_MACHINE_TYPE="n2-standard-8"
export DOWNSTREAM_DISK_GB=200        # normalized betas/qc across all batches + conda env + combine step headroom
export DOWNSTREAM_INSTANCE_NAME="ppmi-downstream-$(date +%Y%m%d-%H%M%S)"
