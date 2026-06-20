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

# ── Load methylation matrix ───────────────────────────────────────────────────
cat(sprintf("[%s] Loading methylation...\n", Sys.time()))
meth_df  <- fread(opt$meth_file, data.table=FALSE)
probe_ids <- meth_df[[1]]
meth_mat  <- as.matrix(meth_df[,-1])
rownames(meth_mat) <- probe_ids
cat(sprintf("Methylation: %d probes x %d samples\n",
    nrow(meth_mat), ncol(meth_mat)))

# ── Load ID map and align samples ─────────────────────────────────────────────
# Methylation columns are Pegid; genotype columns are GWAS_ID
# Use id_map to align
id_map <- read.csv(opt$id_map, stringsAsFactors=FALSE)

# Find samples present in both
geno_gwas_ids <- colnames(geno_mat)
meth_pegids   <- colnames(meth_mat)

aligned <- id_map %>%
  filter(Pegid %in% meth_pegids, GWAS_ID %in% geno_gwas_ids)

cat(sprintf("Samples aligned (meth + geno): %d\n", nrow(aligned)))

if (nrow(aligned)==0) {
  stop("No samples align between methylation and genotype matrices")
}

# Subset and reorder both matrices to aligned samples
meth_aligned <- meth_mat[, aligned$Pegid,  drop=FALSE]
geno_aligned <- geno_mat[, aligned$GWAS_ID, drop=FALSE]

# Rename genotype columns to Pegid for consistent MatrixEQTL input
colnames(geno_aligned) <- aligned$Pegid

# ── Load covariates ───────────────────────────────────────────────────────────
covar_df  <- read.csv(opt$covar_file, stringsAsFactors=FALSE,
                      row.names=1, check.names=FALSE)
covar_mat <- as.matrix(covar_df)
# Strip R's automatic X prefix from column names starting with digits
colnames(covar_mat) <- sub("^X", "", colnames(covar_mat))

# Diagnostic: confirm column names and alignment
cat(sprintf("Covariate matrix: %d vars x %d samples\n",
    nrow(covar_mat), ncol(covar_mat)))
cat(sprintf("First 3 covar col names: %s\n",
    paste(head(colnames(covar_mat), 3), collapse=", ")))
cat(sprintf("First 3 aligned Pegids:  %s\n",
    paste(head(aligned$Pegid, 3), collapse=", ")))
cat(sprintf("First 3 meth col names:  %s\n",
    paste(head(colnames(meth_aligned), 3), collapse=", ")))

# Three-way intersection: meth x geno x covar
pegids_in_covar <- intersect(aligned$Pegid, colnames(covar_mat))
pegids_missing  <- setdiff(aligned$Pegid, colnames(covar_mat))

cat(sprintf("Pegids in covar:  %d\n", length(pegids_in_covar)))
cat(sprintf("Pegids missing:   %d\n", length(pegids_missing)))

if (length(pegids_missing) > 0) {
  cat(sprintf("WARNING: %d Pegids missing from covar — first 5:\n",
              length(pegids_missing)))
  cat(paste(head(pegids_missing, 5), collapse=", "), "\n")
  aligned      <- aligned %>% filter(Pegid %in% pegids_in_covar)
  meth_aligned <- meth_aligned[, aligned$Pegid, drop=FALSE]
  geno_aligned <- geno_aligned[, aligned$Pegid, drop=FALSE]
}

if (nrow(aligned) == 0) {
  stop("No samples remain after three-way alignment — cannot proceed.")
}

# Final covariate subset — guaranteed to succeed
covar_aligned <- covar_mat[, aligned$Pegid, drop=FALSE]

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
