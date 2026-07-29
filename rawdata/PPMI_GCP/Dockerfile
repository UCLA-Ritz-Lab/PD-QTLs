# SeSAMe normalization container for Google Batch tasks.
# Base: Bioconductor image already has R + BiocManager preconfigured.
FROM bioconductor/bioconductor_docker:RELEASE_3_19

# --- gcloud CLI (for gsutil/gcloud storage sync inside the task) ---
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl apt-transport-https ca-certificates gnupg \
    && echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
        | tee /etc/apt/sources.list.d/google-cloud-sdk.list \
    && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
        | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - \
    && apt-get update && apt-get install -y google-cloud-cli \
    && rm -rf /var/lib/apt/lists/*

# --- R packages required by process_batch.R, plus EpiDISH (used later
# downstream, harmless to include here too so the same image can be reused) ---
RUN R -e "BiocManager::install(c('sesame', 'sesameData', 'GenomicRanges', 'EpiDISH'), update = FALSE, ask = FALSE)"
RUN R -e "install.packages(c('tidyverse', 'data.table', 'optparse'), repos = 'https://cloud.r-project.org')"

# Pre-cache SeSAMe reference data at build time so tasks don't each hit the
# network on first use (this is the step that caused download failures
# mid-run when done per-task previously). Note: sesameDataCache() no longer
# accepts a bare platform keyword like "EPIC" -- it wants explicit
# ExperimentHub titles, or sesameDataCacheAll() to grab everything. We use
# the latter here since this is a one-time build-time cost, not a per-task
# cost, and it avoids guessing wrong about which specific titles
# process_batch.R needs internally (e.g. EPIC.address, EPIC.1.SigDF,
# EPIC.addressSpecies, EPIC.probeInfo) and having a task fail mid-fanout on
# a missing title.
RUN R -e "library(sesame); library(sesameData); sesameDataCacheAll()"

WORKDIR /pipeline
COPY scripts/ /pipeline/scripts/
RUN chmod +x /pipeline/scripts/entrypoint.sh

ENTRYPOINT ["/pipeline/scripts/entrypoint.sh"]
