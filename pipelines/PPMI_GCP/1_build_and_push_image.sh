#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source config.sh

gcloud auth configure-docker "${REGION}-docker.pkg.dev" --quiet

# Cloud Build avoids needing Docker installed locally and matches the
# pattern already used for the Pioneer Metals Cloud Run deploys.
gcloud builds submit --tag "$IMAGE_URI" .

echo "Image pushed: $IMAGE_URI"
