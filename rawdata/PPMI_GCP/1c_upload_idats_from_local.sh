#!/usr/bin/env bash
# Uploads local IDATs directly to the bucket. Safe to Ctrl-C and rerun --
# `gcloud storage rsync` compares source/dest and only transfers what's
# missing or changed, and transfers files in parallel by default (no need
# for the old gsutil -m flag).
#
# Usage: ./1c_upload_idats_from_local.sh /path/to/Downloads/idats
set -uo pipefail   # deliberately not -e: the retry loop below handles failures
cd "$(dirname "$0")"
source config.sh

LOCAL_IDAT_DIR="${1:?Usage: $0 /path/to/Downloads/idats}"

echo "=== Uploading ${LOCAL_IDAT_DIR} -> gs://${BUCKET}/${RAW_PREFIX}/ ==="
echo "This is resumable — safe to interrupt (Ctrl-C, connection drop, closed laptop) and rerun."

until gcloud storage rsync "$LOCAL_IDAT_DIR" "gs://${BUCKET}/${RAW_PREFIX}" --recursive; do
  echo "rsync interrupted or failed — retrying in 30s (Ctrl-C to stop entirely)..."
  sleep 30
done

echo "=== Upload complete ==="
