#!/usr/bin/env Rscript
# run_mqtl.R
# Runs MatrixEQTL cis+trans mQTL analysis for one dataset.
# Called per dataset by the run_mqtl Snakemake rule.

suppressPackageStartupMessages({
  library(MatrixEQTL)
  library(data.table)
  library(tidyverse)
  library(optparse)
})

option_list <- list(
  make_option("--dataset",     type="character"),
  make_option("--geno_prefix", type="character"),
  make_option("--meth_file",   type="character"),
  make_option("--covar_file",  type="character"),
  make_option("--id_map",      type="character"),
  make_option("--snp_pos",     type="character"),
  make_option("--cpg_pos",     type="character"),
  make_option("--out_cis",     type="character"),
  make_option("--out_trans",   type="character"),
  make_option("--cis_dist",    type="double",  default=1e6),
  make_option("--pval_cis",    type="double",  default=1e-3),
  make_option("--pval_trans",  type="double",  default=1e-5),
  make_option("--threads",     type="integer", default=4),
  make_option("--probe_list",  type="character", default=NULL,
             help="Optional restricted probe list for intersection analysis"),
  make_option("--log",         type="character", default="run_mqtl.log")
)
opt <- parse_args(OptionParser(option_list=option_list))

log_con <- file(opt$log, open="wt")
sink(log_con, type="message")
on.exit({ sink(type="message"); close(log_con) })

cat(sprintf("[%s] Running mQTL: %s\n", Sys.time(), opt$dataset))

# ── Load genotype dosage matrix ───────────────────────────────────────────────
# PLINK2 exported additive dosage (.raw file)
raw_file <- paste0(opt$geno_prefix, ".raw")
cat(sprintf("[%s] Loading genotypes from %s...\n", Sys.time(), raw_file))
geno_raw <- fread(raw_file, data.table=FALSE)
snp_cols  <- colnames(geno_raw)[-(1:6)]
geno_mat  <- t(as.matrix(geno_raw[, snp_cols]))
colnames(geno_mat) <- geno_raw$IID
# Strip PLINK allele suffix from SNP names (e.g. rs123_A -> rs123)
rownames(geno_mat) <- sub("_[ACGT]$", "", rownames(geno_mat))
cat(sprintf("Genotype matrix: %d SNPs x %d samples\n",
    nrow(geno_mat), ncol(geno_mat)))
cat(sprintf("First 3 raw genotype IIDs: %s\n",
    paste(head(colnames(geno_mat), 3), collapse=", ")))

# ── Load methylation matrix ───────────────────────────────────────────────────
cat(sprintf("[%s] Loading methylation...\n", Sys.time()))
# Load methylation matrix — read header explicitly to avoid fread
# mangling digit-starting column names into V1,V2,V3...
meth_con    <- gzcon(file(opt$meth_file, "rb"))
meth_header <- strsplit(readLines(meth_con, n=1), ",")[[1]]
close(meth_con)

# meth_header[1] = "probe_id", meth_header[2:end] = sample PATNOs
sample_cols <- meth_header[-1]

cat(sprintf("Header sample cols (first 3): %s\n",
    paste(head(sample_cols, 3), collapse=", ")))

# Load full matrix — fread header=TRUE so it reads its own header
# then we overwrite with the correctly parsed one
meth_df <- fread(opt$meth_file, data.table=FALSE, header=TRUE)

# Rows: probes, cols: probe_id + samples
# Use sample_cols (correctly parsed) as column names for sample columns
probe_ids <- meth_df[[1]]
meth_mat  <- as.matrix(meth_df[, -1, drop=FALSE])
rownames(meth_mat) <- probe_ids
colnames(meth_mat) <- sample_cols   # assign correctly parsed names

cat(sprintf("Methylation: %d probes x %d samples\n",
    nrow(meth_mat), ncol(meth_mat)))
cat(sprintf("First 3 meth col names: %s\n",
    paste(head(colnames(meth_mat), 3), collapse=", ")))

# Restrict to intersection probe list if provided
if (!is.null(opt$probe_list) && file.exists(opt$probe_list)) {
  keep_probes <- readLines(opt$probe_list)
  n_before    <- nrow(meth_mat)
  meth_mat    <- meth_mat[rownames(meth_mat) %in% keep_probes, , drop=FALSE]
  cat(sprintf("Probes after intersection filter: %d / %d\n",
              nrow(meth_mat), n_before))
}

# ── Load ID map and align samples ─────────────────────────────────────────────
# Methylation columns are Pegid; genotype columns are GWAS_ID
# Use id_map to align
id_map <- read.csv(opt$id_map, stringsAsFactors=FALSE)

# Find samples present in both
geno_gwas_ids <- colnames(geno_mat)
meth_pegids   <- colnames(meth_mat)

# PPMI id_map uses PATNO (methylation) and plink_iid (genotype)
# The raw file IIDs are in PPMISI{PATNO}.variant2 format — match via id_map
geno_gwas_ids <- colnames(geno_mat)
cat(sprintf("First 3 geno IIDs: %s\n",
    paste(head(geno_gwas_ids, 3), collapse=", ")))
cat(sprintf("First 3 id_map plink_iids: %s\n",
    paste(head(id_map$plink_iid, 3), collapse=", ")))

aligned <- id_map %>%
  filter(PATNO %in% meth_pegids, plink_iid %in% geno_gwas_ids)

cat(sprintf("Samples aligned (meth + geno): %d\n", nrow(aligned)))

if (nrow(aligned)==0) {
  stop("No samples align between methylation and genotype matrices")
}

# Diagnose any remaining mismatches before subsetting
patno_in_meth   <- intersect(aligned$PATNO, colnames(meth_mat))
patno_missing   <- setdiff(aligned$PATNO, colnames(meth_mat))
cat(sprintf("aligned$PATNO in meth_mat: %d\n", length(patno_in_meth)))
cat(sprintf("aligned$PATNO NOT in meth_mat: %d\n", length(patno_missing)))
if (length(patno_missing) > 0) {
  cat(sprintf("First missing: '%s' vs meth col '%s'\n",
              patno_missing[1], colnames(meth_mat)[1]))
  cat(sprintf("nchar missing: %d, nchar meth: %d\n",
              nchar(patno_missing[1]), nchar(colnames(meth_mat)[1])))
  # Only keep aligned rows where PATNO is in meth_mat
  aligned <- aligned %>% filter(PATNO %in% colnames(meth_mat))
  cat(sprintf("Aligned after meth filter: %d\n", nrow(aligned)))
}

# Ensure all IDs are plain character with no factors, NAs or whitespace
aligned$PATNO     <- trimws(as.character(aligned$PATNO))
aligned$plink_iid <- trimws(as.character(aligned$plink_iid))
colnames(meth_mat) <- trimws(colnames(meth_mat))
colnames(geno_mat) <- trimws(colnames(geno_mat))

# Use match() for safe indexed subsetting instead of name-based subscript
meth_idx <- match(aligned$PATNO, colnames(meth_mat))
geno_idx <- match(aligned$plink_iid, colnames(geno_mat))

cat(sprintf("meth_idx NAs: %d\n", sum(is.na(meth_idx))))
cat(sprintf("geno_idx NAs: %d\n", sum(is.na(geno_idx))))

# Remove any unmatched samples
valid <- !is.na(meth_idx) & !is.na(geno_idx)
aligned   <- aligned[valid, ]
meth_idx  <- meth_idx[valid]
geno_idx  <- geno_idx[valid]

cat(sprintf("Samples after index matching: %d\n", nrow(aligned)))

meth_aligned <- meth_mat[, meth_idx, drop=FALSE]
geno_aligned <- geno_mat[, geno_idx, drop=FALSE]
colnames(meth_aligned) <- aligned$PATNO
colnames(geno_aligned) <- aligned$PATNO

# Rename genotype columns to PATNO for consistent MatrixEQTL input
colnames(geno_aligned) <- aligned$PATNO

# ── Load covariates ───────────────────────────────────────────────────────────
covar_df  <- read.csv(opt$covar_file, stringsAsFactors=FALSE,
                      row.names=1, check.names=FALSE)
covar_mat <- as.matrix(covar_df)
# Strip R's automatic X prefix from column names starting with digits
colnames(covar_mat) <- sub("^X", "", colnames(covar_mat))

# Diagnostic: confirm column names and alignment
cat(sprintf("Covariate matrix: %d vars x %d samples\n",
    nrow(covar_mat), ncol(covar_mat)))
cat(sprintf("First 3 covar col names:    %s\n",
    paste(head(colnames(covar_mat), 3), collapse=", ")))
cat(sprintf("First 3 aligned$PATNO:      %s\n",
    paste(head(aligned$PATNO, 3), collapse=", ")))
cat(sprintf("First 3 meth_aligned cols:  %s\n",
    paste(head(colnames(meth_aligned), 3), collapse=", ")))
cat(sprintf("nchar covar col[1]: %d, nchar PATNO[1]: %d\n",
    nchar(colnames(covar_mat)[1]), nchar(aligned$PATNO[1])))

# Three-way intersection: meth x geno x covar
patnos_in_covar <- intersect(aligned$PATNO, colnames(covar_mat))
patnos_missing  <- setdiff(aligned$PATNO, colnames(covar_mat))

cat(sprintf("PATNOs in covar:  %d\n", length(patnos_in_covar)))
cat(sprintf("PATNOs missing:   %d\n", length(patnos_missing)))

if (length(patnos_missing) > 0) {
  cat(sprintf("WARNING: %d PATNOs missing from covar — first 5:\n",
              length(patnos_missing)))
  cat(paste(head(patnos_missing, 5), collapse=", "), "\n")
  aligned      <- aligned %>% filter(PATNO %in% patnos_in_covar)
  meth_aligned <- meth_aligned[, aligned$PATNO, drop=FALSE]
  geno_aligned <- geno_aligned[, aligned$PATNO, drop=FALSE]
}

if (nrow(aligned) == 0) {
  stop("No samples remain after three-way alignment — cannot proceed.")
}

# Final covariate subset
cat(sprintf("Pre-covar subset: aligned n=%d, covar cols=%d\n",
    nrow(aligned), ncol(covar_mat)))
cat(sprintf("PATNOs in covar: %d\n",
    sum(aligned$PATNO %in% colnames(covar_mat))))
covar_aligned <- covar_mat[, aligned$PATNO, drop=FALSE]

# Ensure numeric storage (coercion from data.frame can produce character)
storage.mode(covar_aligned) <- "double"

cat(sprintf("Final aligned samples: %d\n", ncol(covar_aligned)))
cat(sprintf("Covariates: %d vars x %d samples\n",
    nrow(covar_aligned), ncol(covar_aligned)))

# ── Load SNP and CpG positions ────────────────────────────────────────────────
snp_pos <- fread(opt$snp_pos, data.table=FALSE)
cpg_pos <- fread(opt$cpg_pos, data.table=FALSE)

# Filter to SNPs/CpGs present in our matrices
snp_pos <- snp_pos %>% filter(snp_id %in% rownames(geno_aligned))
cpg_pos <- cpg_pos %>% filter(cpg_id %in% rownames(meth_aligned))

cat(sprintf("SNP positions: %d\n", nrow(snp_pos)))
cat(sprintf("CpG positions: %d\n", nrow(cpg_pos)))

# ── Build MatrixEQTL SlicedData objects ───────────────────────────────────────
snp_data   <- SlicedData$new(); snp_data$CreateFromMatrix(geno_aligned)
meth_data  <- SlicedData$new(); meth_data$CreateFromMatrix(meth_aligned)
covar_data <- SlicedData$new(); covar_data$CreateFromMatrix(covar_aligned)

# ── Run MatrixEQTL ────────────────────────────────────────────────────────────
cat(sprintf("[%s] Running MatrixEQTL...\n", Sys.time()))

# MatrixEQTL requires:
#   snpspos: data.frame with columns snpid, chr, pos
#   genepos: data.frame with columns geneid, chr, s1, s2
# CpGs are single-base so s1 == s2
snpspos_df <- snp_pos %>%
  select(snpid=snp_id, chr, pos) %>%
  mutate(chr=as.character(chr))

genepos_df <- cpg_pos %>%
  select(geneid=cpg_id, chr, s1=pos) %>%
  mutate(s2=s1, chr=as.character(chr))

# Filter to SNPs/CpGs actually present in our matrices
snpspos_df <- snpspos_df %>% filter(snpid %in% rownames(geno_aligned))
genepos_df <- genepos_df %>% filter(geneid %in% rownames(meth_aligned))

cat(sprintf("SNP positions for MatrixEQTL: %d\n", nrow(snpspos_df)))
cat(sprintf("CpG positions for MatrixEQTL: %d\n", nrow(genepos_df)))

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

cat(sprintf("[%s] Complete.\n", Sys.time()))
cat(sprintf("Cis mQTLs  (p < %.0e): %d\n", opt$pval_cis,  me$cis$neqtls))
cat(sprintf("Trans mQTLs (p < %.0e): %d\n", opt$pval_trans, me$trans$neqtls))
