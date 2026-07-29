#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source config.sh

echo "=== Building manifest of Sentrix batches under gs://${BUCKET}/${RAW_PREFIX}/ ==="
gcloud storage ls "gs://${BUCKET}/${RAW_PREFIX}/" \
  | grep -E '/$' \
  | sed -E "s#gs://${BUCKET}/${RAW_PREFIX}/##; s#/\$##" \
  | sort > /tmp/manifest.txt

N=$(wc -l < /tmp/manifest.txt)
if [[ "$N" -eq 0 ]]; then
  echo "No Sentrix batch directories found under gs://${BUCKET}/${RAW_PREFIX}/ — upload IDATs first." >&2
  exit 1
fi
echo "Found ${N} Sentrix batches:"
cat /tmp/manifest.txt

gcloud storage cp /tmp/manifest.txt "gs://${BUCKET}/${RAW_PREFIX}/manifest.txt"

PARALLELISM=$(( N < NORMALIZE_MAX_PARALLEL ? N : NORMALIZE_MAX_PARALLEL ))
JOB_NAME="ppmi-normalize-$(date +%Y%m%d-%H%M%S)"

cat > /tmp/batch_job.json <<EOF
{
  "taskGroups": [
    {
      "taskSpec": {
        "runnables": [
          {
            "container": {
              "imageUri": "${IMAGE_URI}"
            }
          }
        ],
        "environment": {
          "variables": {
            "BUCKET": "${BUCKET}",
            "RAW_PREFIX": "${RAW_PREFIX}",
            "NORM_PREFIX": "${NORM_PREFIX}",
            "BETAS_PREFIX": "${BETAS_PREFIX}",
            "ARRAY_TYPE": "${ARRAY_TYPE}"
          }
        },
        "computeResource": {
          "cpuMilli": ${NORMALIZE_CPU_MILLI},
          "memoryMib": ${NORMALIZE_MEM_MIB}
        },
        "maxRetryCount": 2
      },
      "taskCount": ${N},
      "parallelism": ${PARALLELISM}
    }
  ],
  "allocationPolicy": {
    "instances": [
      {
        "policy": {
          "provisioningModel": "SPOT",
          "machineType": "e2-standard-2"
        }
      }
    ],
    "network": {
      "networkInterfaces": [
        {
          "network": "projects/${PROJECT_ID}/global/networks/default",
          "subnetwork": "projects/${PROJECT_ID}/regions/${REGION}/subnetworks/default",
          "noExternalIpAddress": true
        }
      ]
    },
    "serviceAccount": {
      "email": "${BATCH_SA}"
    },
    "location": {
      "allowedLocations": ["zones/${ZONE}"]
    }
  },
  "logsPolicy": {
    "destination": "CLOUD_LOGGING"
  }
}
EOF

echo "=== Submitting ${JOB_NAME}: ${N} tasks, ${PARALLELISM} concurrent Spot VMs ==="
gcloud batch jobs submit "$JOB_NAME" \
  --location="$REGION" \
  --config=/tmp/batch_job.json

echo "Track progress with:"
echo "  gcloud batch jobs describe ${JOB_NAME} --location=${REGION}"
echo "  gcloud logging read 'resource.type=batch.googleapis.com/Job AND resource.labels.job_id=${JOB_NAME}' --limit=50"
