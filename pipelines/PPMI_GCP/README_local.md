# PPMI EPIC v2 Methylation Pipeline: IDATs → Beta Matrix

Processes raw IDAT files from the Illumina Infinium EPIC v2 array into a
clean, normalized probe × sample beta matrix ready for longitudinal mQTL
analysis.

## Pipeline Overview

```
per-batch IDAT directories (Downloads/idats/{SENTRIX_ID}/)
    │
    ├─ Step 1: validate_sample_sheet   check IDATs present, columns valid
    ├─ Step 2: normalize_batch         SeSAMe NOOB + pOOBAH per batch
    ├─ Step 3: convert_rds_to_csv      RDS → CSV.gz for Python
    │
    └─ Step 4: combine_betas           intersect probes, merge all batches
                │
                ├─ Step 5: deconvolve_celltypes   EpiDISH blood proportions
                └─ Step 6: qc_report              summary across all batches
```

## Directory Structure Required

```
Downloads/idats/
├── 203286230019/
│   ├── sample_sheet.csv
│   ├── 203286230019_R01C01_Red.idat
│   ├── 203286230019_R01C01_Grn.idat
│   └── ...
├── 203286230020/
│   └── ...
```

## Sample Sheet Format (per batch)

```
Sample_Name,Sentrix_ID,Sentrix_Position,Sample_Group
3156,203286230019,R05C01,Healthy Control_BL
3157,203286230019,R06C01,Parkinson_V04
```

- `Sample_Name` — subject ID (used as SubjectID in beta matrix columns)
- `Sentrix_ID` — matches the directory name
- `Sentrix_Position` — R{row}C{col} format
- `Sample_Group` — diagnosis + visit code separated by underscore
  - Visit codes: BL, V04, V06, V08, V10, V12

## Setup

### 1. Install Snakemake

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

### 2. Pre-cache SeSAMe Reference Data

SeSAMe downloads array manifests on first use. Run this once before the
pipeline to avoid download failures mid-run:

```bash
conda run -n sesame_env Rscript - << 'EOF'
library(sesame)
library(SeSAMeData)
sesameDataCache("EPICv2")
EOF
```

This downloads ~500 MB of reference data to your R cache directory.

### 3. Edit config.yaml

Key parameters to review:

```yaml
idat_dir: "../../../Downloads/PPMI_GCP/idats"    # path to your IDAT directories
out_dir:  "../../../results/PPMI_GCP"            # output directory

sesame:
  pval_threshold: 0.05         # pOOBAH masking threshold
  min_sample_detection_rate: 0.85  # drop sample if < 85% probes detected
  min_probe_detection_rate: 0.90   # drop probe if < 90% samples detected
```

### 4. Dry Run

```bash
conda activate snakemake
snakemake --use-conda --cores 1 --dry-run
```

Confirm the number of batches matches your IDAT directory count.

### 5. Run the Pipeline

```bash
# Recommended: limit cores to control RAM (2 R threads per batch)
snakemake --use-conda --cores 8
```

With `threads_r: 2` per batch and `--cores 8`, up to 4 batches run in
parallel. Each batch uses ~2–4 GB RAM, so 4 parallel batches = ~16 GB peak,
well within your 32 GB limit.

---

## Output Files

```
results/
├── validated/          per-batch IDAT validation sentinels
├── betas/              per-batch beta matrices (.rds and .csv.gz)
├── qc/                 per-batch QC metrics (.qc.csv)
├── sample_sheets/      per-batch parsed sample sheets with visit info
├── logs/               per-step log files
└── final/
    ├── beta_matrix.csv.gz          final probe × sample matrix  ← main output
    ├── beta_matrix.parquet         same matrix in Parquet format (faster I/O)
    ├── cell_type_proportions.csv   EpiDISH blood proportions per sample
    ├── qc_combined.csv             all-sample QC metrics
    ├── retained_probes.txt         list of probes passing all filters
    └── qc_report.txt               human-readable QC summary
```

## Connecting to the mQTL RMarkdown

In `dynamic_mqtl_workflow.Rmd`, set:

```r
BETA_FILE   <- "results/final/beta_matrix.csv.gz"
SAMPLE_SHEET <- "results/final/qc_combined.csv"
```

The beta matrix columns are named `{SubjectID}_T{visit_num}` matching the
convention expected by the mQTL RMarkdown (e.g. `3156_T1` for subject 3156
at baseline, `3156_T2` for V04, etc.).

## RAM Usage Guide

| Parallel batches | Approx. peak RAM |
|------------------|-----------------|
| 1                | ~3 GB           |
| 2                | ~6 GB           |
| 4                | ~12 GB          |
| 8                | ~24 GB          |

Stay at 4 parallel batches (`--cores 8` with `threads_r: 2`) to keep
comfortably within 32 GB.

## Probe Filters Applied

| Filter                  | Rationale                                      |
|-------------------------|------------------------------------------------|
| pOOBAH p > 0.05         | Poor hybridization signal — unreliable beta    |
| Sex chromosome probes   | Confounded by subject sex composition          |
| Cross-reactive probes   | Map to multiple genomic locations              |
| SNP-confounded probes   | Genotype effects inflate methylation variance  |
| > 10% NA across samples | Unreliable across cohort                       |

## PPMI Visit Code → Visit Number Mapping

| PPMI Code | Visit Num | Description          |
|-----------|-----------|----------------------|
| BL        | 1         | Baseline             |
| V04       | 2         | Visit 4 (6 months)   |
| V06       | 3         | Visit 6 (12 months)  |
| V08       | 4         | Visit 8 (18 months)  |
| V10       | 5         | Visit 10 (24 months) |
| V12       | 6         | Visit 12 (30 months) |
