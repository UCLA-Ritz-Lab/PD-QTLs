#!/usr/bin/env bash
# Runs the sequential, non-parallel part of the pipeline (convert_rds_to_csv
# -> combine_betas -> deconvolve_celltypes -> qc_report) on a single Spot VM,
# then uploads the compact final outputs and deletes itself.
#
# One-time prerequisite: your Snakemake pipeline repo (Snakefile, envs/,
# scripts/) needs to exist at gs://$BUCKET/pipeline/. Do this once with:
#   gsutil -m rsync -r ./ppmi_methylation_pipeline gs://$BUCKET/pipeline/
set -euo pipefail
cd "$(dirname "$0")"
source config.sh

cat > /tmp/startup.sh <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec > >(tee /var/log/startup-script.log) 2>&1

echo "=== Installing Miniforge + Snakemake ==="
# Miniforge instead of Miniconda: ships with conda-forge as its only default
# channel, so it never touches Anaconda's pkgs/main or pkgs/r channels --
# those now require an explicit non-interactive Terms of Service acceptance
# (conda tos accept ...) before use, which would otherwise fail here. We
# only ever install from conda-forge/bioconda anyway, so there's no reason
# to pull in Anaconda's distribution or its ToS gate at all.
curl -sL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -o /tmp/miniforge.sh
bash /tmp/miniforge.sh -b -p /opt/conda
source /opt/conda/etc/profile.d/conda.sh
conda create -y -n snakemake -c conda-forge -c bioconda snakemake py-spy
conda activate snakemake

echo "=== Syncing pipeline repo and normalized batch results ==="
mkdir -p /opt/pipeline
gcloud storage cp -r "gs://${BUCKET}/pipeline/*" /opt/pipeline/
mkdir -p /opt/pipeline/results/normalized /opt/pipeline/results/betas
gcloud storage cp -r "gs://${BUCKET}/${NORM_PREFIX}/*" /opt/pipeline/results/normalized/
# Pull the CSVs the fan-out step already converted inline (see
# entrypoint.sh) -- touched after the RDS files so Snakemake's mtime check
# sees convert_rds_to_csv's output as newer than its input and skips
# re-running it, going straight to combine_betas.
#gcloud storage cp -r "gs://${BUCKET}/${BETAS_PREFIX}/*" /opt/pipeline/results/betas/
#touch /opt/pipeline/results/betas/*.csv.gz

echo "=== Running downstream Snakemake rules ==="
cd /opt/pipeline
snakemake --cores $(nproc) --use-conda \
  results/final/beta_matrix.csv.gz \
  results/final/cell_type_proportions.csv \
  results/final/qc_report.txt \
  results/final/retained_probes.txt

echo "=== Uploading compact final outputs ==="
gcloud storage cp -r \
  results/final/beta_matrix.csv.gz \
  results/final/cell_type_proportions.csv \
  results/final/qc_report.txt \
  results/final/retained_probes.txt \
  "gs://${BUCKET}/${FINAL_PREFIX}/"

echo "=== Done — self-deleting instance ==="
INSTANCE_NAME=\$(curl -s -H "Metadata-Flavor: Google" "http://metadata.google.internal/computeMetadata/v1/instance/name")
ZONE_FULL=\$(curl -s -H "Metadata-Flavor: Google" "http://metadata.google.internal/computeMetadata/v1/instance/zone")
ZONE_SHORT=\$(basename "\$ZONE_FULL")
gcloud compute instances delete "\$INSTANCE_NAME" --zone="\$ZONE_SHORT" --quiet
EOF

echo "=== Launching Spot VM ${DOWNSTREAM_INSTANCE_NAME} ==="
# Unlike the Batch normalize tasks, this VM's startup script needs general
# internet access (curl to repo.anaconda.com for Miniconda), not just
# Google APIs -- so it keeps its external IP rather than using Private
# Google Access. Options to fully avoid external IPs everywhere would be a
# Cloud NAT gateway, but that's unnecessary complexity for a single one-off
# VM.
gcloud compute instances create "$DOWNSTREAM_INSTANCE_NAME" \
  --zone="$ZONE" \
  --machine-type="$DOWNSTREAM_MACHINE_TYPE" \
  --boot-disk-size="${DOWNSTREAM_DISK_GB}GB" \
  --boot-disk-type=pd-balanced \
  --provisioning-model=SPOT \
  --instance-termination-action=STOP \
  --image-family=debian-12 \
  --image-project=debian-cloud \
  --scopes=cloud-platform \
  --service-account="$VM_SA" \
  --metadata-from-file=startup-script=/tmp/startup.sh

echo "VM launched. It will run unattended and self-delete when finished."
echo "Watch progress with:"
echo "  gcloud compute instances tail-serial-port-output ${DOWNSTREAM_INSTANCE_NAME} --zone=${ZONE}"
echo "Once you see the instance is gone (gcloud compute instances list), pull results with:"
echo "  gcloud storage cp -r gs://${BUCKET}/${FINAL_PREFIX} ./local_results/"
