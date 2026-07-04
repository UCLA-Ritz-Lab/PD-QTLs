#!/usr/bin/env Rscript
# prepare_datasets.R
# Prepares PPMI baseline mQTL input files:
#   - Extracts T1 (baseline) methylation columns
#   - Matches subject IDs across methylation, genotypes, phenotype, covariates
#   - Stratifies by case/control status
#   - Outputs per-dataset methylation, covariate, and genotype keep files

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--beta_matrix",     type="character"),
  make_option("--cell_type_props", type="character"),
  make_option("--plink_prefix",    type="character"),
  make_option("--phenotype_file",  type="character"),
  make_option("--out_dir",         type="character"),
  make_option("--meth_visit",      type="integer",   default=1L),
  make_option("--plink_id_prefix", type="character", default="PPMISI"),
  make_option("--plink_id_suffix", type="character", default=".variant2"),
  make_option("--cohort_case",     type="integer",   default=1L),
  make_option("--cohort_control",  type="integer",   default=2L),
  make_option("--min_sd",          type="double",    default=0.02),
  make_option("--max_na_frac",     type="double",    default=0.05),
  make_option("--log",             type="character", default="prepare_datasets.log")
)
opt <- parse_args(OptionParser(option_list=option_list))

dir.create(dirname(opt$log), recursive=TRUE, showWarnings=FALSE)
log_con <- file(opt$log, open="wt")
sink(log_con, type="message")
sink(log_con, type="output", append=TRUE, split=TRUE)
on.exit({ sink(type="output"); sink(type="message"); close(log_con) })

log_section <- function(t) {
  cat("\n", strrep("=",60), "\n", t, "\n", strrep("=",60), "\n", sep="")
}

# ── Step 1: Load phenotype file ───────────────────────────────────────────────
log_section("Step 1: Load phenotype file")

pheno <- read.csv(opt$phenotype_file, stringsAsFactors=FALSE) %>%
  select(PATNO, COHORT, COHORT_DEFINITION, ENROLL_AGE) %>%
  mutate(PATNO = as.character(PATNO))

cat(sprintf("Phenotype file: %d subjects\n", nrow(pheno)))
cat("COHORT distribution:\n")
print(table(pheno$COHORT, pheno$COHORT_DEFINITION))

cases    <- pheno %>% filter(COHORT == opt$cohort_case)
controls <- pheno %>% filter(COHORT == opt$cohort_control)
cat(sprintf("Cases (COHORT=%d):    %d\n", opt$cohort_case,    nrow(cases)))
cat(sprintf("Controls (COHORT=%d): %d\n", opt$cohort_control, nrow(controls)))

# ── Step 2: Load PLINK2 sample IDs ───────────────────────────────────────────
log_section("Step 2: Load PLINK2 sample IDs")

psam <- read.table(
  paste0(opt$plink_prefix, ".psam"),
  header=TRUE, comment.char="", stringsAsFactors=FALSE
)
colnames(psam)[1] <- "IID"   # handle #IID header

cat(sprintf("PLINK2 samples (raw): %d\n", nrow(psam)))

# Keep only .variant2 entries — these contain actual genotype data
# .variant entries exist in psam but have no genotype information
psam_primary <- psam %>%
  filter(grepl(paste0(opt$plink_id_suffix, "$"), IID, fixed=FALSE)) %>%
  mutate(
    PATNO = str_remove(IID, paste0("^", opt$plink_id_prefix)),
    PATNO = str_remove(PATNO, paste0(opt$plink_id_suffix, "$"))
  )

cat(sprintf("PLINK2 samples (primary .variant2 only): %d\n", nrow(psam_primary)))
cat(sprintf("First 5 IID -> PATNO mappings:\n"))
print(head(psam_primary %>% select(IID, PATNO), 5))

# ── Step 3: Load baseline methylation ────────────────────────────────────────
log_section("Step 3: Load baseline methylation")

cat("Loading beta matrix header to identify baseline columns...\n")

# Read header only first to find T1 columns
header_line <- readLines(gzcon(file(opt$beta_matrix, "rb")), n=1)
all_cols    <- strsplit(header_line, ",")[[1]]

# Find baseline columns: {PATNO}_T{visit}
visit_suffix  <- paste0("_T", opt$meth_visit)
baseline_cols <- all_cols[grepl(paste0(visit_suffix, "$"), all_cols)]
baseline_cols <- baseline_cols[baseline_cols != ""]  # remove empty first col

cat(sprintf("Total columns in beta matrix: %d\n", length(all_cols)))
cat(sprintf("Baseline (%s) columns: %d\n", visit_suffix, length(baseline_cols)))

# Extract PATNO from column names
baseline_patnos <- str_remove(baseline_cols, paste0(visit_suffix, "$"))
cat(sprintf("First 5 baseline PATNOs: %s\n",
    paste(head(baseline_patnos, 5), collapse=", ")))

# ── Step 4: Three-way sample intersection ─────────────────────────────────────
log_section("Step 4: Sample intersection")

# Build master sample table
master <- data.frame(
  PATNO      = baseline_patnos,
  meth_col   = baseline_cols,
  stringsAsFactors = FALSE
) %>%
  inner_join(psam_primary %>% select(PATNO, plink_iid=IID), by="PATNO") %>%
  inner_join(pheno %>% select(PATNO, COHORT, COHORT_DEFINITION, ENROLL_AGE),
             by="PATNO")

cat(sprintf("Methylation baseline samples:    %d\n", length(baseline_patnos)))
cat(sprintf("PLINK2 primary samples:          %d\n", nrow(psam_primary)))
cat(sprintf("Phenotype subjects:              %d\n", nrow(pheno)))
cat(sprintf("Three-way intersection:          %d\n", nrow(master)))
cat(sprintf("  Cases:    %d\n", sum(master$COHORT == opt$cohort_case)))
cat(sprintf("  Controls: %d\n", sum(master$COHORT == opt$cohort_control)))

# ── Step 5: Load full baseline beta matrix ────────────────────────────────────
log_section("Step 5: Load baseline beta matrix")

cat("Loading baseline columns from beta matrix (may take a few minutes)...\n")

# Read only the baseline columns + probe ID column
# Column indices to read (1-based, including first unnamed probe ID column)
col_indices <- which(all_cols %in% c("", baseline_cols))

beta_full <- fread(
  opt$beta_matrix,
  select    = c(1L, which(all_cols %in% baseline_cols)),
  data.table= FALSE
)

# Set probe IDs as rownames
probe_ids  <- beta_full[[1]]
beta_mat   <- as.matrix(beta_full[, -1])
rownames(beta_mat) <- probe_ids

# Rename columns from {PATNO}_T1 to PATNO for clean downstream use
colnames(beta_mat) <- str_remove(colnames(beta_mat), paste0(visit_suffix, "$"))

cat(sprintf("Baseline beta matrix: %d probes x %d samples\n",
    nrow(beta_mat), ncol(beta_mat)))

# ── Step 6: Probe filtering ───────────────────────────────────────────────────
log_section("Step 6: Probe filtering")

probe_sds <- apply(beta_mat, 1, sd,        na.rm=TRUE)
probe_nas <- rowMeans(is.na(beta_mat))

keep_probes <- rownames(beta_mat)[
  probe_sds >= opt$min_sd & probe_nas <= opt$max_na_frac
]
cat(sprintf("Removed (SD < %.2f):        %d\n",
    opt$min_sd,     sum(probe_sds < opt$min_sd)))
cat(sprintf("Removed (NA > %.0f%%):       %d\n",
    100*opt$max_na_frac, sum(probe_nas > opt$max_na_frac & probe_sds >= opt$min_sd)))
cat(sprintf("Probes retained:             %d\n", length(keep_probes)))

beta_filt <- beta_mat[keep_probes, ]

# Save retained probe list
dir.create(file.path(opt$out_dir, "datasets"), recursive=TRUE, showWarnings=FALSE)
writeLines(keep_probes,
           file.path(opt$out_dir, "datasets", "retained_probes.txt"))

# ── Step 7: Load cell-type proportions ───────────────────────────────────────
log_section("Step 7: Load cell-type proportions")

cell_props <- read.csv(opt$cell_type_props, stringsAsFactors=FALSE)

# sample_col_name format is {PATNO}_T{visit} — extract PATNO
cell_props <- cell_props %>%
  mutate(
    PATNO = str_remove(sample_col_name, paste0(visit_suffix, "$"))
  ) %>%
  filter(grepl(paste0(visit_suffix, "$"), sample_col_name))

cat(sprintf("Cell-type proportions at baseline: %d samples\n", nrow(cell_props)))
cat(sprintf("Cell-type columns: %s\n",
    paste(setdiff(colnames(cell_props),
                  c("sample_col_name","PATNO")), collapse=", ")))

# ── Step 8: Build covariate matrix ───────────────────────────────────────────
log_section("Step 8: Build covariate matrix")

# Covariates: age + cell-type proportions
# Drop one cell type to avoid collinearity (Neutro is typically most abundant)
cell_cols <- setdiff(colnames(cell_props),
                     c("sample_col_name", "PATNO", "Neutro", "Eosino"))

covar_data <- master %>%
  select(PATNO, ENROLL_AGE) %>%
  left_join(cell_props %>% select(PATNO, all_of(cell_cols)), by="PATNO")

n_missing_covar <- sum(!master$PATNO %in% covar_data$PATNO[
  complete.cases(covar_data)
])
if (n_missing_covar > 0) {
  cat(sprintf("WARNING: %d samples missing covariates\n", n_missing_covar))
}

cat(sprintf("Covariate columns: %s\n",
    paste(c("ENROLL_AGE", cell_cols), collapse=", ")))

# ── Step 9: Save per-dataset outputs ─────────────────────────────────────────
log_section("Step 9: Save datasets")

save_dataset <- function(cohort_val, dataset_name, label) {
  cat(sprintf("\n[%s] %s\n", dataset_name, label))

  # Samples for this stratum — join master with complete covariate rows
  # Use suffix to handle ENROLL_AGE appearing in both master and covar_data
  stratum <- master %>%
    filter(COHORT == cohort_val) %>%
    select(PATNO, plink_iid, COHORT, COHORT_DEFINITION, ENROLL_AGE) %>%
    inner_join(
      covar_data %>%
        filter(complete.cases(.)) %>%
        select(PATNO, all_of(cell_cols)),   # only cell cols from covar_data
      by = "PATNO"
    )

  cat(sprintf("  Samples: %d\n", nrow(stratum)))
  cat(sprintf("  Columns: %s\n", paste(colnames(stratum), collapse=", ")))

  if (nrow(stratum) == 0) {
    cat("  WARNING: no samples — skipping\n")
    return(invisible(NULL))
  }

  out <- file.path(opt$out_dir, "datasets", dataset_name)
  dir.create(out, recursive=TRUE, showWarnings=FALSE)

  # Methylation: probes x samples (PATNO columns)
  meth_sub <- beta_filt[, stratum$PATNO, drop=FALSE]
  mf <- file.path(out, "methylation.csv.gz")
  fwrite(cbind(probe_id=rownames(meth_sub), as.data.frame(meth_sub)), mf)
  cat(sprintf("  Methylation: %d probes x %d samples -> %s\n",
              nrow(meth_sub), ncol(meth_sub), basename(mf)))

  # Covariate matrix: covariates x samples (PATNO columns)
  covar_cols_out <- c("ENROLL_AGE", cell_cols)
  covar_sub_df   <- stratum %>% select(PATNO, all_of(covar_cols_out))

  # Verify all covariate columns exist
  missing_cv <- setdiff(covar_cols_out, colnames(covar_sub_df))
  if (length(missing_cv) > 0) {
    stop(sprintf("Missing covariate columns in stratum: %s",
                 paste(missing_cv, collapse=", ")))
  }

  covar_sub <- t(as.matrix(covar_sub_df[, covar_cols_out]))
  colnames(covar_sub) <- stratum$PATNO
  rownames(covar_sub) <- covar_cols_out
  storage.mode(covar_sub) <- "double"

  cf <- file.path(out, "covariates.csv")
  write.csv(cbind(covariate=rownames(covar_sub), as.data.frame(covar_sub)),
            cf, row.names=FALSE, quote=FALSE)
  cat(sprintf("  Covariates: %d vars x %d samples -> %s\n",
              nrow(covar_sub), ncol(covar_sub), basename(cf)))

  # PLINK2 --keep file: single column of PLINK2 IIDs
  gf <- file.path(out, "geno_keep.txt")
  writeLines(stratum$plink_iid, gf)
  cat(sprintf("  Geno keep: %d samples -> %s\n", nrow(stratum), basename(gf)))

  # ID map: PATNO <-> plink_iid
  write.csv(stratum %>% select(PATNO, plink_iid),
            file.path(out, "id_map.csv"),
            row.names=FALSE, quote=FALSE)

  # Sample metadata
  write.csv(stratum, file.path(out, "sample_metadata.csv"),
            row.names=FALSE, quote=FALSE)

  invisible(list(n_samples=nrow(stratum), n_probes=nrow(meth_sub)))
}

result_cases    <- save_dataset(opt$cohort_case,    "ppmi_cases",    "PPMI PD Cases")
result_controls <- save_dataset(opt$cohort_control, "ppmi_controls", "PPMI Controls")

log_section("Summary")
cat(sprintf("%-20s  %8s  %8s\n", "Dataset", "Samples", "Probes"))
cat(strrep("-", 42), "\n")
for (nm in c("ppmi_cases","ppmi_controls")) {
  r <- get(paste0("result_", sub("ppmi_","",nm)))
  if (!is.null(r))
    cat(sprintf("%-20s  %8d  %8d\n", nm, r$n_samples, r$n_probes))
}
cat("\nOutput:", file.path(opt$out_dir, "datasets"), "\n")
