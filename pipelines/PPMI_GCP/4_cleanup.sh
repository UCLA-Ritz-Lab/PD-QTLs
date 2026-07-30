#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source config.sh

echo "=== Confirm final outputs are downloaded locally before proceeding ==="
read -p "Have you already pulled gs://${BUCKET}/${FINAL_PREFIX}/ to your local machine? [y/N] " CONFIRM
if [[ "$CONFIRM" != "y" && "$CONFIRM" != "Y" ]]; then
  echo "Aborting — run: gcloud storage cp -r gs://${BUCKET}/${FINAL_PREFIX} ./local_results/  first."
  exit 1
fi

echo "=== Deleting raw IDATs and normalized intermediates from the bucket ==="
gcloud storage rm -r "gs://${BUCKET}/${RAW_PREFIX}" || true
gcloud storage rm -r "gs://${BUCKET}/${NORM_PREFIX}" || true

echo "=== Checking for any VMs still running ==="
gcloud compute instances list --filter="zone:(${ZONE})"

echo "=== Checking for any Batch jobs still active ==="
gcloud batch jobs list --location="$REGION" --filter="status.state=RUNNING"

echo "Raw/intermediate data removed. Bucket now only holds ${FINAL_PREFIX}/ and the container image in Artifact Registry."
echo "To tear down everything (irreversible), run:"
echo "  gcloud projects delete ${PROJECT_ID}"
