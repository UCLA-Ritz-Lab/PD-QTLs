#!/usr/bin/env Rscript
# run_mqtl.R (PPMI)
# Runs MatrixEQTL cis+trans mQTL analysis for one dataset.
# Based on the working PPMI run_mqtl_old.R with the following changes:
#   - --geno_prefix/.raw replaced with --geno_tsv (SNP-major bgzipped TSV)
#   - CreateFromMatrix replaced with file-based SlicedData$LoadFile
#   - rm()/gc() at every safe point to minimise peak RAM
#   - noFDRsaveMemory=TRUE to avoid accumulating all p-values in RAM

suppressPackageStartupMessages({
  library(MatrixEQTL)
  library(data.table)
  library(tidyverse)
  library(optparse)
})

option_list <- list(
  make_option("--dataset",     type="character"),
  make_option("--geno_tsv",    type="character",
              help="SNP-major bgzipped TSV from export_genotypes rule"),
  make_option("--meth_file",   type="character"),
  make_option("--covar_file",  type="character"),
  make_option("--id_map",      type="character"),
  make_option("--snp_pos",     type="character"),
  make_option("--cpg_pos",     type="character"),
  make_option("--out_cis",     type="character"),
  make_option("--out_trans",   type="character"),
  make_option("--cis_dist",    type="double",  default=1e6),
  make_option("--pval_cis",    type="double",  default=1e-5),
  make_option("--pval_trans",  type="double",  default=1e-8),
  make_option("--threads",     type="integer", default=1),
  make_option("--log",         type="character", default="run_mqtl.log")
)
opt <- parse_args(OptionParser(option_list=option_list))

log_con <- file(opt$log, open="wt")
sink(log_con, type="message")
sink(log_con, type="output", append=TRUE, split=TRUE)
on.exit({ sink(type="output"); sink(type="message"); close(log_con) })

cat(sprintf("[%s] Running mQTL: %s\n", Sys.time(), opt$dataset))

# ── Read genotype sample IDs from header only ─────────────────────────────────
# geno_tsv is SNP-major (rows=SNPs, cols=samples), bgzipped.
# Only the header is read here for alignment; the full matrix is loaded
# later via SlicedData$LoadFile to avoid the fread+t() segfault on large files.
cat(sprintf("[%s] Reading genotype header: %s\n", Sys.time(), opt$geno_tsv))
geno_header  <- fread(cmd=paste("zcat", opt$geno_tsv, "| head -1"),
                      header=FALSE, data.table=FALSE)
geno_samples <- trimws(as.character(geno_header[1, -1]))  # drop snp_id column
cat(sprintf("Genotype file: %d samples in header\n", length(geno_samples)))
cat(sprintf("First 3 geno IIDs: %s\n",
    paste(head(geno_samples, 3), collapse=", ")))
rm(geno_header); gc()

# ── Load methylation matrix ───────────────────────────────────────────────────
cat(sprintf("[%s] Loading methylation...\n", Sys.time()))

# Read header explicitly to avoid fread mangling digit-starting column names
meth_con    <- gzcon(file(opt$meth_file, "rb"))
meth_header <- strsplit(readLines(meth_con, n=1), ",")[[1]]
close(meth_con)
sample_cols <- meth_header[-1]

cat(sprintf("Header sample cols (first 3): %s\n",
    paste(head(sample_cols, 3), collapse=", ")))

meth_df  <- fread(opt$meth_file, data.table=FALSE, header=TRUE)
probe_ids <- meth_df[[1]]
meth_mat  <- as.matrix(meth_df[, -1, drop=FALSE])
rownames(meth_mat) <- probe_ids
colnames(meth_mat) <- sample_cols
rm(meth_df); gc()

cat(sprintf("Methylation: %d probes x %d samples\n",
    nrow(meth_mat), ncol(meth_mat)))
cat(sprintf("First 3 meth col names: %s\n",
    paste(head(colnames(meth_mat), 3), collapse=", ")))

# ── Load ID map and align samples ─────────────────────────────────────────────
# id_map has columns PATNO (methylation ID) and plink_iid (genotype ID)
id_map <- read.csv(opt$id_map, stringsAsFactors=FALSE)

meth_pegids <- trimws(colnames(meth_mat))
colnames(meth_mat) <- meth_pegids

cat(sprintf("First 3 id_map plink_iids: %s\n",
    paste(head(id_map$plink_iid, 3), collapse=", ")))

aligned <- id_map %>%
  filter(PATNO %in% meth_pegids, plink_iid %in% geno_samples)

cat(sprintf("Samples aligned (meth + geno): %d\n", nrow(aligned)))

if (nrow(aligned) == 0) {
  stop("No samples align between methylation and genotype matrices")
}

# Diagnose any remaining mismatches
patno_in_meth <- intersect(aligned$PATNO, colnames(meth_mat))
patno_missing <- setdiff(aligned$PATNO, colnames(meth_mat))
cat(sprintf("aligned$PATNO in meth_mat: %d\n", length(patno_in_meth)))
cat(sprintf("aligned$PATNO NOT in meth_mat: %d\n", length(patno_missing)))
if (length(patno_missing) > 0) {
  cat(sprintf("First missing: '%s' vs meth col '%s'\n",
              patno_missing[1], colnames(meth_mat)[1]))
  aligned <- aligned %>% filter(PATNO %in% colnames(meth_mat))
  cat(sprintf("Aligned after meth filter: %d\n", nrow(aligned)))
}

# Ensure clean character IDs
aligned$PATNO     <- trimws(as.character(aligned$PATNO))
aligned$plink_iid <- trimws(as.character(aligned$plink_iid))

# Use match() for safe indexed subsetting
meth_idx <- match(aligned$PATNO,     colnames(meth_mat))
geno_idx <- match(aligned$plink_iid, geno_samples)

cat(sprintf("meth_idx NAs: %d\n", sum(is.na(meth_idx))))
cat(sprintf("geno_idx NAs: %d\n", sum(is.na(geno_idx))))

valid    <- !is.na(meth_idx) & !is.na(geno_idx)
aligned  <- aligned[valid, ]
meth_idx <- meth_idx[valid]
geno_idx <- geno_idx[valid]

cat(sprintf("Samples after index matching: %d\n", nrow(aligned)))

meth_aligned <- meth_mat[, meth_idx, drop=FALSE]
colnames(meth_aligned) <- aligned$PATNO
rm(meth_mat); gc()

# ── Load covariates ───────────────────────────────────────────────────────────
covar_df  <- read.csv(opt$covar_file, stringsAsFactors=FALSE,
                      row.names=1, check.names=FALSE)
covar_mat <- as.matrix(covar_df)
colnames(covar_mat) <- sub("^X", "", colnames(covar_mat))
rm(covar_df); gc()

cat(sprintf("Covariate matrix: %d vars x %d samples\n",
    nrow(covar_mat), ncol(covar_mat)))
cat(sprintf("First 3 covar col names: %s\n",
    paste(head(colnames(covar_mat), 3), collapse=", ")))
cat(sprintf("First 3 aligned$PATNO:   %s\n",
    paste(head(aligned$PATNO, 3), collapse=", ")))

# Three-way intersection: meth x geno x covar
patnos_in_covar <- intersect(aligned$PATNO, colnames(covar_mat))
patnos_missing  <- setdiff(aligned$PATNO,   colnames(covar_mat))

cat(sprintf("PATNOs in covar:  %d\n", length(patnos_in_covar)))
cat(sprintf("PATNOs missing:   %d\n", length(patnos_missing)))

if (length(patnos_missing) > 0) {
  cat(sprintf("WARNING: %d PATNOs missing from covar — first 5:\n",
              length(patnos_missing)))
  cat(paste(head(patnos_missing, 5), collapse=", "), "\n")
  aligned      <- aligned %>% filter(PATNO %in% patnos_in_covar)
  meth_aligned <- meth_aligned[, aligned$PATNO, drop=FALSE]
}

if (nrow(aligned) == 0) {
  stop("No samples remain after three-way alignment — cannot proceed.")
}

covar_aligned <- covar_mat[, aligned$PATNO, drop=FALSE]
storage.mode(covar_aligned) <- "double"
rm(covar_mat); gc()

cat(sprintf("Final aligned samples: %d\n", ncol(covar_aligned)))
cat(sprintf("Covariates: %d vars x %d samples\n",
    nrow(covar_aligned), ncol(covar_aligned)))

# ── Load SNP and CpG positions ────────────────────────────────────────────────
snp_pos <- fread(opt$snp_pos, data.table=FALSE)
cpg_pos <- fread(opt$cpg_pos, data.table=FALSE)

# ── Write methylation and covariates to temp files ────────────────────────────
# Genotypes: pre-filter geno_tsv to MAF-filtered SNPs before SlicedData loading.
# The full geno_tsv has ~13.7M SNPs; loading all into SlicedData causes OOM.
# snp_pos contains only SNPs passing MAF >= 0.05 (~3-4M), reducing the file
# substantially before MatrixEQTL further restricts to cis-window SNPs.
tmp_meth   <- tempfile(fileext=".txt")
tmp_covar  <- tempfile(fileext=".txt")
tmp_geno   <- tempfile(fileext=".tsv.gz")
tmp_snpids <- tempfile(fileext=".txt")
on.exit(unlink(c(tmp_meth, tmp_covar, tmp_geno, tmp_snpids)), add=TRUE)

write.table(meth_aligned,  tmp_meth,  sep="\t", quote=FALSE)
write.table(covar_aligned, tmp_covar, sep="\t", quote=FALSE)
# Note: meth_aligned freed after column count assertion below; covar_aligned freed now
rm(covar_aligned); gc()

writeLines(snp_pos$snp_id, tmp_snpids)
cat(sprintf("[%s] Pre-filtering geno_tsv to %d MAF-filtered SNPs...\n",
            Sys.time(), nrow(snp_pos)))

filter_cmd <- paste0(
  "{ zcat ", opt$geno_tsv, " | head -1; ",
  "zcat ", opt$geno_tsv, " | tail -n +2 | grep -Fw -f ", tmp_snpids, "; }",
  " | bgzip -c > ", tmp_geno
)
ret <- system(filter_cmd)
if (ret != 0) stop("Failed to pre-filter geno_tsv")

n_filtered <- as.integer(system(
  paste("zcat", tmp_geno, "| tail -n +2 | wc -l"), intern=TRUE))
cat(sprintf("[%s] Filtered geno_tsv: %d SNPs retained\n", Sys.time(), n_filtered))
gc()

# ── Subset and reorder geno_tsv columns to aligned samples ───────────────────
# Critical: MatrixEQTL matches matrices by column POSITION not by ID name.
# The geno_tsv may have more samples than meth_aligned if --mind removed
# samples from geno during subset_genotypes but not from methylation.
# We must extract exactly aligned$plink_iid columns in the correct order.
cat(sprintf("[%s] Subsetting geno_tsv to %d aligned samples...\n",
            Sys.time(), nrow(aligned)))

geno_aligned_cols <- match(aligned$plink_iid, geno_samples)
n_missing <- sum(is.na(geno_aligned_cols))
if (n_missing > 0) {
  cat(sprintf("WARNING: %d plink_iids not found in geno_tsv — dropping\n",
              n_missing))
  keep              <- !is.na(geno_aligned_cols)
  aligned           <- aligned[keep, ]
  geno_aligned_cols <- geno_aligned_cols[keep]
}

# awk 1-based col indices: col 1 = snp_id, sample cols start at 2
awk_cols    <- c(1, geno_aligned_cols + 1)
awk_col_str <- paste(paste0("$", awk_cols), collapse=", ")

tmp_geno_aligned <- tempfile(fileext=".tsv.gz")
on.exit(unlink(tmp_geno_aligned), add=TRUE)

col_cmd <- paste0(
  "zcat ", tmp_geno, " | ",
  "awk 'BEGIN{OFS=\"\\t\"}{print ", awk_col_str, "}' | ",
  "bgzip -c > ", tmp_geno_aligned
)
if (system(col_cmd) != 0) stop("Failed to subset geno_tsv columns")

n_geno_cols <- as.integer(system(
  paste("zcat", tmp_geno_aligned, "| head -1 | awk '{print NF-1}'"),
  intern=TRUE))
cat(sprintf("[%s] Geno_tsv: %d SNPs x %d samples\n",
            Sys.time(), n_filtered, n_geno_cols))
cat(sprintf("Meth aligned:  %d samples\n", ncol(meth_aligned)))
stopifnot(n_geno_cols == ncol(meth_aligned))
rm(meth_aligned); gc()

# ── Build MatrixEQTL SlicedData objects ───────────────────────────────────────
snp_data   <- SlicedData$new()
meth_data  <- SlicedData$new()
covar_data <- SlicedData$new()

snp_data$fileDelimiter      <- "\t"
snp_data$fileOmitCharacters <- "NA"
snp_data$fileSkipRows       <- 1
snp_data$fileSkipColumns    <- 1
snp_data$fileSliceSize      <- 500
snp_data$LoadFile(tmp_geno_aligned)

meth_data$fileDelimiter      <- "\t"
meth_data$fileOmitCharacters <- "NA"
meth_data$fileSkipRows       <- 1
meth_data$fileSkipColumns    <- 1
meth_data$fileSliceSize      <- 1000
meth_data$LoadFile(tmp_meth)

covar_data$fileDelimiter      <- "\t"
covar_data$fileOmitCharacters <- "NA"
covar_data$fileSkipRows       <- 1
covar_data$fileSkipColumns    <- 1
covar_data$fileSliceSize      <- 1e6
covar_data$LoadFile(tmp_covar)

# ── Build position data frames (after SlicedData loaded) ──────────────────────
snpspos_df <- snp_pos %>%
  filter(snp_id %in% snp_data$GetAllRowNames()) %>%
  select(snpid=snp_id, chr, pos) %>%
  mutate(chr=as.character(chr))

genepos_df <- cpg_pos %>%
  filter(cpg_id %in% meth_data$GetAllRowNames()) %>%
  select(geneid=cpg_id, chr, s1=pos) %>%
  mutate(s2=s1, chr=as.character(chr))

cat(sprintf("SNP positions for MatrixEQTL: %d\n", nrow(snpspos_df)))
cat(sprintf("CpG positions for MatrixEQTL: %d\n", nrow(genepos_df)))

cat(sprintf("[%s] Memory before MatrixEQTL: %.1f MB\n",
    Sys.time(), sum(gc()[,2]) * 8 / 1024))

# ── Run MatrixEQTL ────────────────────────────────────────────────────────────
cat(sprintf("[%s] Running MatrixEQTL...\n", Sys.time()))

me <- Matrix_eQTL_main(
  snps                  = snp_data,
  gene                  = meth_data,
  cvrt                  = covar_data,
  output_file_name      = opt$out_trans,
  pvOutputThreshold     = opt$pval_trans,
  useModel              = modelLINEAR,
  errorCovariance       = numeric(),
  verbose               = TRUE,
  output_file_name.cis  = opt$out_cis,
  pvOutputThreshold.cis = opt$pval_cis,
  snpspos               = snpspos_df,
  genepos               = genepos_df,
  cisDist               = opt$cis_dist,
  pvalue.hist           = FALSE,
  min.pv.by.genesnp     = FALSE,
  noFDRsaveMemory       = FALSE
)

rm(snpspos_df, genepos_df, snp_pos, cpg_pos); gc()
cat(sprintf("[%s] Complete.\n", Sys.time()))
cat(sprintf("Cis mQTLs  (p < %.0e): %d\n", opt$pval_cis,  me$cis$neqtls))

# pval_trans = 0 disables the trans pass outright: MatrixEQTL skips the
# computation (rather than merely filtering the output, which is what any
# non-zero threshold does), returns no $trans element, and never creates the
# trans output file. Snakemake declares that file as a rule output, so write a
# header-only stub to keep the rule's contract. Header matches MatrixEQTL's own
# so downstream readers need no special case.
if (is.null(me$trans)) {
  cat(sprintf("Trans mQTLs: analysis disabled (pval_trans=%g)\n", opt$pval_trans))
  writeLines("SNP\tgene\tbeta\tt-stat\tp-value\tFDR", opt$out_trans)
  cat(sprintf("Wrote header-only trans file: %s\n", opt$out_trans))
} else {
  cat(sprintf("Trans mQTLs (p < %.0e): %d\n", opt$pval_trans, me$trans$neqtls))
}
