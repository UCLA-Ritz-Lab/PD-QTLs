#!/usr/bin/env bash
# Provisions a small STANDARD (non-Spot) VM for pulling IDATs from LONI IDA.
#
# Why not Spot here specifically: IDA requires downloads to be started and
# completed from the same IP that logged in, and times out on inactivity.
# A Spot preemption mid-download would invalidate that session and force a
# restart, same as a Colab disconnect would. This step is I/O-bound, not
# compute-bound, so a small standard VM costs very little regardless.
set -euo pipefail
cd "$(dirname "$0")"
source config.sh

export DOWNLOAD_INSTANCE_NAME="ppmi-idat-download-$(date +%Y%m%d-%H%M%S)"
export DOWNLOAD_DISK_GB=300   # buffer for unzip + upload cycles, not the full 172GB at once

gcloud compute instances create "$DOWNLOAD_INSTANCE_NAME" \
  --zone="$ZONE" \
  --machine-type=e2-medium \
  --boot-disk-size="${DOWNLOAD_DISK_GB}GB" \
  --boot-disk-type=pd-standard \
  --image-family=debian-12 \
  --image-project=debian-cloud \
  --scopes=cloud-platform \
  --service-account="$VM_SA"

echo "VM ${DOWNLOAD_INSTANCE_NAME} created in zone ${ZONE}."
echo ""
echo "Next (manual — this part needs your IDA login, so it can't be scripted):"
echo "  1. gcloud compute ssh ${DOWNLOAD_INSTANCE_NAME} --zone=${ZONE}"
echo "  2. sudo apt-get update && sudo apt-get install -y default-jre unzip"
echo "  3. On ida.loni.usc.edu (from your own browser, this part is fine anywhere):"
echo "       - Build your Genetic Data / IDAT collection as usual"
echo "       - Choose 'Advanced Download' -> download the small URL LIST file"
echo "         (not the data itself) and the IDA Downloader jar"
echo "       - scp both to the VM: gcloud compute scp url_list.txt IDA_downloader.jar \\"
echo "           ${DOWNLOAD_INSTANCE_NAME}:~ --zone=${ZONE}"
echo "  4. On the VM, run the IDA Downloader against that URL list — this is what"
echo "     creates a fresh authenticated session FROM THE VM'S IP, satisfying"
echo "     IDA's same-IP requirement (enter your IDA credentials when it prompts):"
echo "       java -jar IDA_downloader.jar   # add its documented flags for the"
echo "                                       # url list + destination dir — check"
echo "                                       # the usage text the jar prints, or"
echo "                                       # the instructions shown next to the"
echo "                                       # download link on the IDA site"
echo "  5. Stage to the bucket in batches to keep local disk usage small, e.g. per"
echo "     zip file: unzip -> gcloud storage cp -r <Sentrix_ID dirs> gs://${BUCKET}/${RAW_PREFIX}/ -> rm -rf locally"
echo "  6. When done: gcloud compute instances delete ${DOWNLOAD_INSTANCE_NAME} --zone=${ZONE}"
