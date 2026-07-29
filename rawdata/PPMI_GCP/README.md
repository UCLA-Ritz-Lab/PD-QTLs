# PPMI IDAT Pipeline on Google Batch

Runs the `normalize_batch` step of the PPMI methylation pipeline (SeSAMe NOOB +
pOOBAH, per Sentrix batch) as fanned-out Spot tasks on Google Batch, then
finishes the sequential/aggregation steps (`convert_rds_to_csv` →
`combine_betas` → `deconvolve_celltypes` → `qc_report`) on a single
short-lived Spot VM. Raw IDATs never touch your local disk; only the final
beta matrix, cell-type proportions, and QC report come back down.

## Order of operations

```
0_setup_gcp.sh              # one-time: project, APIs, bucket, service accounts, lifecycle rule
      │
      ▼
(manual) download IDATs from LONI IDA → gsutil -m cp -r to gs://$BUCKET/raw/idats/
      │
      ▼
1_build_and_push_image.sh   # builds the SeSAMe + gcloud CLI container, pushes to Artifact Registry
      │
      ▼
2_submit_normalize_batch.sh # fans out one Spot task per Sentrix_ID batch on Google Batch
      │
      ▼
3_run_downstream_vm.sh      # single Spot VM: combine/deconvolve/QC, uploads compact results, self-deletes
      │
      ▼
gsutil -m cp -r gs://$BUCKET/results/final ./local_results/
      │
      ▼
4_cleanup.sh                # deletes raw IDATs from the bucket, removes the image, confirms no VMs left running
```

## Config

Every script sources `config.sh` — edit that first (project ID, region,
bucket name). Keeping this in its own GCP project (not `pioneer-metals-501918`)
is recommended so billing/IAM stay cleanly separated from the Pioneer Metals work.

## What `process_batch.R` actually produces

Each Batch task calls your real `process_batch.R` (openSesame `prep="QCDPB"`,
per-sample QC on detection rate, sex-probe removal, high-NA-rate probe
filtering) for one Sentrix batch, writing into a **shared** results root
(namespaced by filename, so many tasks can write into the same prefix without
colliding):

```
gs://$BUCKET/results/normalized/
├── betas/{sentrix_id}.betas.rds
├── qc/{sentrix_id}.qc.csv
├── sample_sheets/{sentrix_id}.parsed.csv
└── logs/{sentrix_id}.log
```

Note the script expects `--idat_dir` to be the *parent* of the Sentrix_ID
folder (it does `file.path(idat_dir, sentrix_id)` internally) — the
entrypoint handles that distinction, so raw IDATs stay laid out in the
bucket exactly as `Downloads/idats/{SENTRIX_ID}/...` was locally.

`process_batch.R` defaults `--array_type` to `EPICv2`, but this PPMI batch is
actually **EPIC v1** — `config.sh` sets `ARRAY_TYPE="EPIC"` and the Dockerfile
pre-caches the EPIC v1 manifest accordingly, overriding the script's own
default via the CLI flag.

## Cost notes

- Spot VMs / Spot Batch tasks: ~60–91% off on-demand.
- Google Batch itself has no service fee — you only pay for the Compute Engine
  resources the tasks consume.
- The lifecycle rule in `0_setup_gcp.sh` auto-deletes anything under `raw/` after
  14 days as a backstop in case `4_cleanup.sh` is never run.
- Egress: only the final compact outputs (`results/final/`) should ever leave
  GCP — expect low tens of MB, well under the free egress tier.
