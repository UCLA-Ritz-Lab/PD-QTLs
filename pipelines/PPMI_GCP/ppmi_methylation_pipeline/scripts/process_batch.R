#!/usr/bin/env Rscript
# process_batch.R
# Per-batch IDAT processing with SeSAMe
# Called by Snakemake normalize_batch rule — one invocation per Sentrix batch
#
# Outputs per batch:
#   {out_dir}/betas/{sentrix_id}.betas.rds     — probe x sample beta matrix
#   {out_dir}/qc/{sentrix_id}.qc.csv           — per-sample QC metrics
#   {out_dir}/sample_sheets/{sentrix_id}.parsed.csv — parsed sample metadata

suppressPackageStartupMessages({
  library(sesame)
  library(sesameData)
  library(GenomicRanges)
  library(tidyverse)
  library(data.table)
  library(optparse)
})

# Restrict R to single-threaded BLAS/LAPACK to avoid pthread_create errors
# when called from Snakemake with limited thread resources
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")
Sys.setenv(NUMEXPR_NUM_THREADS = "1")

# ── Argument parsing ───────────────────────────────────────────────────────────
option_list <- list(
  make_option("--sentrix_id",        type = "character"),
  make_option("--idat_dir",          type = "character"),
  make_option("--out_dir",           type = "character"),
  make_option("--array_type",        type = "character", default = "EPICv2"),
  make_option("--pval_threshold",    type = "double",    default = 0.05),
  make_option("--min_sample_rate",   type = "double",    default = 0.85),
  make_option("--min_probe_rate",    type = "double",    default = 0.90),
  make_option("--remove_sex",        type = "logical",   default = TRUE),
  make_option("--remove_cr",         type = "logical",   default = TRUE),
  make_option("--remove_snp",        type = "logical",   default = TRUE),
  make_option("--normalization",     type = "character", default = "NOOB"),
  make_option("--log",               type = "character", default = "process_batch.log")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Redirect messages and warnings to log file
log_con <- file(opt$log, open = "wt")
sink(log_con, type = "message")
on.exit({ sink(type = "message"); close(log_con) })

cat(sprintf("[%s] Processing batch: %s\n", Sys.time(), opt$sentrix_id))

# ── Define output paths early so they are available for early-exit conditions ──
out_betas <- file.path(opt$out_dir, "betas",
                        paste0(opt$sentrix_id, ".betas.rds"))
out_qc    <- file.path(opt$out_dir, "qc",
                        paste0(opt$sentrix_id, ".qc.csv"))
out_ss    <- file.path(opt$out_dir, "sample_sheets",
                        paste0(opt$sentrix_id, ".parsed.csv"))

dir.create(dirname(out_betas), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_qc),    recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_ss),    recursive = TRUE, showWarnings = FALSE)

# ── Locate IDAT files and sample sheet ────────────────────────────────────────
batch_dir <- file.path(opt$idat_dir, opt$sentrix_id)

sample_sheet_path <- file.path(batch_dir, "sample_sheet.csv")
if (!file.exists(sample_sheet_path)) {
  stop(sprintf("sample_sheet.csv not found in %s", batch_dir))
}

# Parse sample sheet
# Columns: Sample_Name, Sentrix_ID, Sentrix_Position, Sample_Group
ss <- read_csv(sample_sheet_path, show_col_types = FALSE) %>%
  mutate(
    sentrix_id       = as.character(Sentrix_ID),
    sentrix_position = Sentrix_Position,
    # Construct full IDAT prefix path: dir/SENTRIX_ID_POSITION
    idat_prefix      = file.path(batch_dir,
                                  paste0(sentrix_id, "_", sentrix_position)),
    # Parse subject ID and visit from Sample_Name and Sample_Group
    subject_id       = as.character(Sample_Name),
    # Visit code is after the last underscore in Sample_Group
    # e.g. "Healthy Control_BL" -> "BL", "Parkinson_V04" -> "V04"
    visit_code       = str_extract(Sample_Group, "[^_]+$"),
    diagnosis        = str_extract(Sample_Group, "^[^_]+"),
    # Convert PPMI visit codes to numeric order for longitudinal modeling
    # PPMI visit codes mapped to sequential integers preserving temporal order.
    # Even-numbered visits (V04, V06...) are standard 6-month intervals.
    # Odd-numbered visits (V03, V05, V07, V09) appear in some cohorts at
    # intermediate timepoints (~3-month intervals).
    visit_num        = case_when(
      visit_code == "BL"  ~ 1L,
      visit_code == "V03" ~ 2L,
      visit_code == "V04" ~ 3L,
      visit_code == "V05" ~ 4L,
      visit_code == "V06" ~ 5L,
      visit_code == "V07" ~ 6L,
      visit_code == "V08" ~ 7L,
      visit_code == "V09" ~ 8L,
      visit_code == "V10" ~ 9L,
      visit_code == "V12" ~ 10L,
      visit_code == "V14" ~ 11L,
      visit_code == "V16" ~ 12L,
      visit_code == "V18" ~ 13L,
      visit_code == "V20" ~ 14L,
      visit_code == "ST"  ~ 15L,
      visit_code == "SC"  ~ 16L,
      TRUE                ~ NA_integer_
    ),
    # Construct column name for beta matrix: SubjectID_T{visit_num}
    sample_col_name  = paste0(subject_id, "_T", visit_num),
    batch            = opt$sentrix_id
  )

cat(sprintf("[%s] Found %d samples in batch\n", Sys.time(), nrow(ss)))

# Verify all IDAT files exist
missing_idats <- ss %>%
  filter(
    !file.exists(paste0(idat_prefix, "_Red.idat")) |
    !file.exists(paste0(idat_prefix, "_Grn.idat"))
  )
if (nrow(missing_idats) > 0) {
  cat(sprintf("[WARN] Missing IDAT files for %d samples:\n", nrow(missing_idats)))
  cat(paste(missing_idats$idat_prefix, collapse = "\n"), "\n")
}

ss <- ss %>%
  filter(
    file.exists(paste0(idat_prefix, "_Red.idat")),
    file.exists(paste0(idat_prefix, "_Grn.idat"))
  )

cat(sprintf("[%s] %d samples with valid IDAT pairs\n", Sys.time(), nrow(ss)))

# ── Load IDATs via SeSAMe ──────────────────────────────────────────────────────
# openSesame processes each sample sequentially, keeping RAM flat
# regardless of batch size

cat(sprintf("[%s] Loading IDATs with SeSAMe (%s)...\n",
            Sys.time(), opt$array_type))

# Build named vector of IDAT prefixes: names = sample column names
idat_prefixes <- setNames(ss$idat_prefix, ss$sample_col_name)

# openSesame pipeline:
#   - prepSesame: loads raw intensities
#   - NOOB: background correction + dye-bias normalisation
#   - pOOBAH: detection p-value masking (sets failed probes to NA)
#   - getBetas: convert to beta values

sdfs <- lapply(seq_along(idat_prefixes), function(i) {
  sname  <- names(idat_prefixes)[i]
  prefix <- idat_prefixes[i]
  cat(sprintf("[%s]   Loading %s (%d/%d)\n",
              Sys.time(), sname, i, length(idat_prefixes)))

  tryCatch({
    sdf <- readIDATpair(prefix, platform = opt$array_type)
    sdf
  }, error = function(e) {
    cat(sprintf("[ERROR] Failed to load %s: %s\n", sname, e$message))
    # Retry without platform specification as fallback
    tryCatch({
      cat(sprintf("[RETRY] Attempting %s without explicit platform...\n", sname))
      sdf <- readIDATpair(prefix)
      sdf
    }, error = function(e2) {
      cat(sprintf("[ERROR] Retry also failed for %s: %s\n", sname, e2$message))
      NULL
    })
  })
})
names(sdfs) <- names(idat_prefixes)

# Remove failed samples
failed_load <- names(sdfs)[sapply(sdfs, is.null)]
if (length(failed_load) > 0) {
  cat(sprintf("[WARN] %d samples failed to load: %s\n",
              length(failed_load), paste(failed_load, collapse = ", ")))
}
sdfs <- Filter(Negate(is.null), sdfs)

# ── Per-sample QC metrics ──────────────────────────────────────────────────────
# Compute detection rate directly from pOOBAH p-values on the raw SigDF.
# sesameQC_calcStats is not used here as it requires a preprocessed SigDF
# and its internal pOOBAH call fails on raw signal in newer SeSAMe versions.
cat(sprintf("[%s] Computing QC metrics...\n", Sys.time()))

qc_metrics <- lapply(names(sdfs), function(sname) {
  sdf <- sdfs[[sname]]

  # pOOBAH returns a SigDF with pval column when return.pval = FALSE (default)
  # Get p-values by running pOOBAH and extracting the pval slot
  tryCatch({
    pvals <- pOOBAH(sdf, return.pval = TRUE)
    detection_rate <- sum(pvals <= opt$pval_threshold, na.rm = TRUE) /
                      length(pvals)

    # Mean intensity from raw channels as a simple QC metric
    mean_int <- mean(c(sdf$MG, sdf$MR, sdf$UG, sdf$UR), na.rm = TRUE)

    data.frame(
      sample_col_name  = sname,
      detection_rate   = detection_rate,
      mean_intensity   = mean_int,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat(sprintf("[WARN] QC metrics failed for %s: %s\n", sname, e$message))
    data.frame(
      sample_col_name  = sname,
      detection_rate   = NA_real_,
      mean_intensity   = NA_real_,
      stringsAsFactors = FALSE
    )
  })
}) %>% bind_rows()

# Flag samples failing detection rate threshold
qc_metrics <- qc_metrics %>%
  mutate(
    pass_qc = detection_rate >= opt$min_sample_rate,
    fail_reason = case_when(
      detection_rate < opt$min_sample_rate ~
        sprintf("detection_rate %.3f < threshold %.3f",
                detection_rate, opt$min_sample_rate),
      TRUE ~ NA_character_
    )
  )

n_fail <- sum(!qc_metrics$pass_qc)
cat(sprintf("[%s] QC: %d/%d samples pass (%.1f%%)\n",
            Sys.time(),
            sum(qc_metrics$pass_qc), nrow(qc_metrics),
            100 * mean(qc_metrics$pass_qc)))

if (n_fail > 0) {
  cat("[WARN] Failing samples:\n")
  qc_metrics %>%
    filter(!pass_qc) %>%
    select(sample_col_name, detection_rate, fail_reason) %>%
    print()
}

# ── Normalize and extract betas via openSesame ────────────────────────────────
# openSesame with prep="QCDPB" runs the full recommended pipeline:
#   Q = channel inference (corrects red/green swap)
#   C = NOOB background correction
#   D = dye bias correction (type I/II normalisation)
#   P = pOOBAH detection masking (sets failed probes to NA)
#   B = getBetas (convert to beta values)
# qualityMask is applied separately after to mask cross-reactive/SNP probes.
# This is cleaner and more robust than the manual step-by-step approach.
cat(sprintf("[%s] Normalizing and extracting beta values (openSesame)...\n",
            Sys.time()))

# Keep only QC-passing samples
passing_samples <- qc_metrics %>% filter(pass_qc) %>% pull(sample_col_name)
pass_prefixes   <- idat_prefixes[passing_samples]

beta_list <- lapply(seq_along(pass_prefixes), function(i) {
  sname  <- names(pass_prefixes)[i]
  prefix <- pass_prefixes[i]
  cat(sprintf("[%s]   Normalizing %s (%d/%d)\n",
              Sys.time(), sname, i, length(pass_prefixes)))
  tryCatch({
    # Full pipeline in one call
    # BPPARAM = SerialParam() forces single-threaded execution
    # avoiding pthread resource conflicts under Snakemake
    betas <- openSesame(
      prefix,
      platform = opt$array_type,
      prep     = "QCDPB",
      BPPARAM  = BiocParallel::SerialParam()
    )
    betas
  }, error = function(e) {
    cat(sprintf("[ERROR] openSesame failed for %s: %s\n", sname, e$message))
    NULL
  })
})
names(beta_list) <- names(pass_prefixes)
beta_list <- Filter(Negate(is.null), beta_list)

# Guard against empty batch — no samples passed QC
if (length(beta_list) == 0) {
  cat(sprintf("[WARN] No samples passed QC in batch %s — writing empty outputs\n",
              opt$sentrix_id))
  # Write empty placeholder files so Snakemake considers the rule complete
  saveRDS(matrix(nrow = 0, ncol = 0), out_betas)
  write_csv(qc_metrics, out_qc)
  write_csv(ss %>%
              left_join(qc_metrics %>%
                          select(sample_col_name, detection_rate,
                                 pass_qc, fail_reason),
                        by = "sample_col_name"),
            out_ss)
  cat(sprintf("[%s] Batch %s complete (0 samples passed QC).\n",
              Sys.time(), opt$sentrix_id))
  quit(status = 0)
}

# Combine into probe x sample matrix
beta_mat <- do.call(cbind, beta_list)
colnames(beta_mat) <- names(beta_list)
cat(sprintf("[%s] Beta matrix: %d probes x %d samples\n",
            Sys.time(), nrow(beta_mat), ncol(beta_mat)))

# ── Probe-level filtering ──────────────────────────────────────────────────────
# SeSAMe's qualityMask() applies the unified KYCG.EPICv2.Mask object which
# covers cross-reactive, SNP-confounded, and other problematic probes in one
# step. Sex chromosome probes are identified via sesameData_getProbesByRegion().
# These masks are applied to each SigDF before getBetas(), then we also apply
# a post-hoc NA rate filter on the assembled beta matrix.
cat(sprintf("[%s] Applying probe filters...\n", Sys.time()))
n_probes_start <- nrow(beta_mat)

# 1. Remove sex chromosome probes
# sesameData_getProbesByRegion can return NULL for some platform/genome
# combinations. Instead, use the EPIC manifest annotation directly to
# identify probes on chrX and chrY — more reliable across SeSAMe versions.
if (opt$remove_sex) {
  tryCatch({
    # Get manifest with chromosomal coordinates
    manifest <- sesameData_getManifestGRanges(opt$array_type)
    if (!is.null(manifest)) {
      sex_probes <- names(manifest)[
        as.character(GenomicRanges::seqnames(manifest)) %in% c("chrX", "chrY")
      ]
    } else {
      # Fallback: use the address annotation from sesameData
      addr <- sesameDataGet(paste0(opt$array_type, ".address"))
      # Filter by chromosome annotation if available
      if ("CpG_chrm" %in% colnames(addr$Manifest)) {
        sex_probes <- rownames(addr$Manifest)[
          addr$Manifest$CpG_chrm %in% c("chrX", "chrY")
        ]
      } else {
        cat("  [WARN] Could not identify sex probes from manifest — skipping\n")
        sex_probes <- character(0)
      }
    }
    n_before <- nrow(beta_mat)
    beta_mat <- beta_mat[!rownames(beta_mat) %in% sex_probes, ]
    cat(sprintf("  Sex probes removed: %d -> %d probes\n",
                n_before, nrow(beta_mat)))
  }, error = function(e) {
    cat(sprintf("  [WARN] Sex probe removal failed: %s\n", e$message))
    cat("  [WARN] Skipping sex probe removal for this batch.\n")
  })
}

# 2. qualityMask (cross-reactive + SNP probes) is already applied by
# openSesame prep="QCDPB" — no separate step needed here.
# Log how many NAs remain from masking for QC transparency.
if (opt$remove_cr || opt$remove_snp) {
  n_masked <- sum(is.na(beta_mat))
  n_total  <- prod(dim(beta_mat))
  cat(sprintf("  Probes masked by openSesame qualityMask: %d/%d (%.1f%%)\n",
              n_masked, n_total, 100 * n_masked / n_total))
}

# 3. Remove probes with high NA rate across this batch
na_rate  <- rowMeans(is.na(beta_mat))
n_before <- nrow(beta_mat)
beta_mat <- beta_mat[na_rate <= (1 - opt$min_probe_rate), ]
cat(sprintf("  High-NA probes removed: %d -> %d probes\n",
            n_before, nrow(beta_mat)))

cat(sprintf("[%s] Final beta matrix: %d probes x %d samples\n",
            Sys.time(), nrow(beta_mat), ncol(beta_mat)))

# ── Save outputs ───────────────────────────────────────────────────────────────
saveRDS(beta_mat, out_betas)
write_csv(qc_metrics, out_qc)
write_csv(ss %>%
            left_join(qc_metrics %>%
                        select(sample_col_name, detection_rate,
                               pass_qc, fail_reason),
                      by = "sample_col_name"),
          out_ss)

cat(sprintf("[%s] Batch %s complete.\n", Sys.time(), opt$sentrix_id))
cat(sprintf("  Beta matrix:   %s\n", out_betas))
cat(sprintf("  QC metrics:    %s\n", out_qc))
cat(sprintf("  Sample sheet:  %s\n", out_ss))
