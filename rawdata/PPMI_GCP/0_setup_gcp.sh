#!/usr/bin/env bash
# One-time setup: project, APIs, bucket + lifecycle rule, Artifact Registry,
# service accounts and IAM bindings.
set -euo pipefail
cd "$(dirname "$0")"
source config.sh

echo "=== 1. Project ==="
# Comment out if the project already exists.
gcloud projects create "$PROJECT_ID" --name="PPMI IDAT Batch Pipeline" || true
gcloud config set project "$PROJECT_ID"

# Link a billing account (replace with your billing account ID — find it via
# `gcloud billing accounts list`).
# gcloud billing projects link "$PROJECT_ID" --billing-account=XXXXXX-XXXXXX-XXXXXX

echo "=== 2. Enable APIs ==="
gcloud services enable \
  batch.googleapis.com \
  compute.googleapis.com \
  storage.googleapis.com \
  artifactregistry.googleapis.com \
  cloudbuild.googleapis.com

echo "=== 3. Bucket ==="
gcloud storage buckets create "gs://${BUCKET}" \
  --location="$REGION" \
  --uniform-bucket-level-access

# Lifecycle rule: auto-delete anything under raw/ after 14 days, as a backstop
# for whenever 4_cleanup.sh doesn't get run manually.
cat > /tmp/lifecycle.json <<EOF
{
  "rule": [
    {
      "action": {"type": "Delete"},
      "condition": {"age": 14, "matchesPrefix": ["${RAW_PREFIX}/"]}
    }
  ]
}
EOF
gcloud storage buckets update "gs://${BUCKET}" --lifecycle-file=/tmp/lifecycle.json

echo "=== 4. Artifact Registry ==="
gcloud artifacts repositories create "$AR_REPO" \
  --repository-format=docker \
  --location="$REGION" \
  --description="PPMI pipeline container images" || true

echo "=== 5. Service accounts ==="
gcloud iam service-accounts create batch-runner \
  --display-name="PPMI Batch task runner" || true
gcloud iam service-accounts create downstream-vm \
  --display-name="PPMI downstream aggregation VM" || true

# Batch tasks: read raw IDATs, write normalized results.
gcloud storage buckets add-iam-policy-binding "gs://${BUCKET}" \
  --member="serviceAccount:${BATCH_SA}" \
  --role="roles/storage.objectAdmin"

# Batch tasks also need to pull the container image from Artifact Registry
# and write their own logs to Cloud Logging -- without these, tasks fail
# silently with no logs anywhere to explain why (image pull failure with no
# logging permission to report it).
gcloud projects add-iam-policy-binding "$PROJECT_ID" \
  --member="serviceAccount:${BATCH_SA}" \
  --role="roles/artifactregistry.reader"
gcloud projects add-iam-policy-binding "$PROJECT_ID" \
  --member="serviceAccount:${BATCH_SA}" \
  --role="roles/logging.logWriter"

# Downstream VM: read normalized results, write final compact outputs, and
# permission to delete itself when done (self-shutdown pattern in step 3).
gcloud storage buckets add-iam-policy-binding "gs://${BUCKET}" \
  --member="serviceAccount:${VM_SA}" \
  --role="roles/storage.objectAdmin"
gcloud projects add-iam-policy-binding "$PROJECT_ID" \
  --member="serviceAccount:${VM_SA}" \
  --role="roles/compute.instanceAdmin.v1" \
  --condition=None
gcloud projects add-iam-policy-binding "$PROJECT_ID" \
  --member="serviceAccount:${VM_SA}" \
  --role="roles/logging.logWriter"

echo "=== 6. Cloud Build default service account ==="
# Recent projects run Cloud Build as the Compute Engine default service
# account (PROJECT_NUMBER-compute@developer.gserviceaccount.com) rather than
# the old dedicated Cloud Build SA, and it has no storage/logging/registry
# access by default -- without these grants, `gcloud builds submit` fails
# reading back its own uploaded source tarball, and would fail again pushing
# the finished image to Artifact Registry.
PROJECT_NUMBER=$(gcloud projects describe "$PROJECT_ID" --format='value(projectNumber)')
COMPUTE_SA="${PROJECT_NUMBER}-compute@developer.gserviceaccount.com"

gcloud projects add-iam-policy-binding "$PROJECT_ID" \
  --member="serviceAccount:${COMPUTE_SA}" \
  --role="roles/storage.objectViewer"
gcloud projects add-iam-policy-binding "$PROJECT_ID" \
  --member="serviceAccount:${COMPUTE_SA}" \
  --role="roles/logging.logWriter"
gcloud projects add-iam-policy-binding "$PROJECT_ID" \
  --member="serviceAccount:${COMPUTE_SA}" \
  --role="roles/artifactregistry.writer"

echo "=== 7. Private Google Access (lets Batch/VM tasks skip external IPs) ==="
# Batch task VMs and the downstream VM only need to reach Google APIs
# (Cloud Storage, Artifact Registry, Cloud Logging) -- not the public
# internet. Without this, each task VM grabs an external IP and you'll hit
# the default IN_USE_ADDRESSES quota (8 per region) as soon as parallelism
# exceeds it. This also avoids the small per-task IAM/network surface of
# giving every task VM a public IP.
gcloud compute networks subnets update default \
  --region="$REGION" \
  --enable-private-ip-google-access

echo "Setup complete. Next: download IDATs from LONI IDA and upload with:"
echo "  gsutil -m cp -r /path/to/local/idats/* gs://${BUCKET}/${RAW_PREFIX}/"
