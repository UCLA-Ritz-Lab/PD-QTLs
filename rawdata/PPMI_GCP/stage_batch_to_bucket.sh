#!/usr/bin/env bash
# Run this ON the download VM (not locally) after `source config.sh` there too.
# Usage: ./stage_batch_to_bucket.sh path/to/downloaded_package.zip
#
# Keeps local disk usage small by processing one IDA zip at a time: unzip,
# push any Sentrix_ID directories it contains up to the bucket, then delete
# the local copies before moving to the next zip.
set -euo pipefail
source config.sh

ZIP_PATH="${1:?Usage: $0 path/to/package.zip}"
WORKDIR=$(mktemp -d)

echo "=== Unzipping $ZIP_PATH ==="
unzip -q "$ZIP_PATH" -d "$WORKDIR"

echo "=== Uploading Sentrix batch directories to gs://${BUCKET}/${RAW_PREFIX}/ ==="
# Adjust this find depth/pattern if your zip's internal layout differs —
# it assumes {SENTRIX_ID}/ directories containing the .idat pairs and
# sample_sheet.csv, matching the local Downloads/idats/{SENTRIX_ID}/ layout.
find "$WORKDIR" -mindepth 1 -maxdepth 3 -type d -regextype posix-extended \
     -regex '.*/[0-9]{9,12}$' | while read -r sentrix_dir; do
  sentrix_id=$(basename "$sentrix_dir")
  echo "  -> ${sentrix_id}"
  gcloud storage cp -r "$sentrix_dir" "gs://${BUCKET}/${RAW_PREFIX}/${sentrix_id}"
done

echo "=== Cleaning up local copies ==="
rm -rf "$WORKDIR"
rm -f "$ZIP_PATH"

echo "Done with $(basename "$ZIP_PATH")."
