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
  make_option("--geno_tsv",    type="character"),
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
  make_option("--gzip_output", action="store_true", default=FALSE),
  make_option("--threads",     type="integer", default=4),
  make_option("--log",         type="character", default="run_mqtl.log")
)
opt <- parse_args(OptionParser(option_list=option_list))

log_con <- file(opt$log, open="wt")
sink(log_con, type="message")
on.exit({ sink(type="message"); close(log_con) })

cat(sprintf("[%s] Running mQTL: %s\n", Sys.time(), opt$dataset))

# ── Read genotype sample IDs from header only ────────────────────────────────
# geno_tsv is SNP-major (rows=SNPs, cols=samples), bgzipped.
# Only the header is read here for sample alignment; the full matrix is
# loaded later via SlicedData$LoadFile to avoid holding it all in RAM
# (the .raw equivalent was 5.2GB and caused a segfault on fread + t()).
cat(sprintf("[%s] Reading genotype header from %s\n", Sys.time(), opt$geno_tsv))
geno_header  <- fread(cmd=paste("zcat", opt$geno_tsv, "| head -1"),
                      header=FALSE, data.table=FALSE)
geno_samples <- as.character(geno_header[1, -1])  # drop snp_id column
cat(sprintf("Genotype file: %d samples in header\n", length(geno_samples)))
rm(geno_header); gc()

# ── Load methylation matrix ───────────────────────────────────────────────────
cat(sprintf("[%s] Loading methylation...\n", Sys.time()))
meth_df  <- fread(opt$meth_file, data.table=FALSE)
probe_ids <- meth_df[[1]]
meth_mat  <- as.matrix(meth_df[,-1])
rownames(meth_mat) <- probe_ids
rm(meth_df); gc()
cat(sprintf("Methylation: %d probes x %d samples\n",
    nrow(meth_mat), ncol(meth_mat)))

# ── Load ID map and align samples ─────────────────────────────────────────────
# Methylation columns are Pegid; geno_samples are GWAS_IDs from tsv header.
# Use id_map to align -- geno matrix itself is not loaded into RAM here.
id_map <- read.csv(opt$id_map, stringsAsFactors=FALSE)

meth_pegids <- colnames(meth_mat)

aligned <- id_map %>%
  filter(Pegid %in% meth_pegids, GWAS_ID %in% geno_samples)

cat(sprintf("Samples aligned (meth + geno): %d\n", nrow(aligned)))

if (nrow(aligned) == 0) {
  stop("No samples align between methylation and genotype matrices")
}

meth_aligned <- meth_mat[, aligned$Pegid, drop=FALSE]
rm(meth_mat); gc()
# Note: geno matrix alignment is enforced when writing the SlicedData-ready
# file below -- geno_samples retained for covariate intersection

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
pegids_in_covar <- Reduce(intersect, list(aligned$Pegid, geno_samples, colnames(covar_mat)))
pegids_missing  <- setdiff(aligned$Pegid, colnames(covar_mat))

cat(sprintf("Pegids in covar:  %d\n", length(pegids_in_covar)))
cat(sprintf("Pegids missing:   %d\n", length(pegids_missing)))

if (length(pegids_missing) > 0) {
  cat(sprintf("WARNING: %d Pegids missing from covar — first 5:\n",
              length(pegids_missing)))
  cat(paste(head(pegids_missing, 5), collapse=", "), "\n")
  aligned      <- aligned %>% filter(Pegid %in% pegids_in_covar)
  meth_aligned <- meth_aligned[, aligned$Pegid, drop=FALSE]
  geno_samples <- aligned$GWAS_ID  # keep geno_samples in sync with aligned
}

if (nrow(aligned) == 0) {
  stop("No samples remain after three-way alignment — cannot proceed.")
}

# Final covariate subset — guaranteed to succeed
covar_aligned <- covar_mat[, aligned$Pegid, drop=FALSE]

# Ensure numeric storage (coercion from data.frame can produce character)
storage.mode(covar_aligned) <- "double"
rm(covar_df, covar_mat); gc()

cat(sprintf("Final aligned samples: %d\n", ncol(covar_aligned)))
cat(sprintf("Covariates: %d vars x %d samples\n",
    nrow(covar_aligned), ncol(covar_aligned)))

# ── Load SNP and CpG positions ────────────────────────────────────────────────
snp_pos <- fread(opt$snp_pos, data.table=FALSE)
cpg_pos <- fread(opt$cpg_pos, data.table=FALSE)

# ── Write methylation and covariates to temp files ────────────────────────────
# Genotypes stay on disk (geno_tsv) and are loaded directly via SlicedData.
tmp_meth  <- tempfile(fileext = ".txt")
tmp_covar <- tempfile(fileext = ".txt")
on.exit(unlink(c(tmp_meth, tmp_covar)), add = TRUE)

write.table(meth_aligned,  tmp_meth,  sep = "\t", quote = FALSE)
write.table(covar_aligned, tmp_covar, sep = "\t", quote = FALSE)
rm(meth_aligned, covar_aligned); gc()

# ── Build MatrixEQTL SlicedData objects ───────────────────────────────────────
snp_data   <- SlicedData$new()
meth_data  <- SlicedData$new()
covar_data <- SlicedData$new()

# Load genotypes directly from bgzipped SNP-major file -- never fully in RAM
snp_data$fileDelimiter      <- "\t"
snp_data$fileOmitCharacters <- "NA"
snp_data$fileSkipRows       <- 1
snp_data$fileSkipColumns    <- 1
snp_data$fileSliceSize      <- 500
snp_data$LoadFile(opt$geno_tsv)

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

# ── Build position data frames (after SlicedData loaded -- rownames available)─
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

# ── Memory state before MatrixEQTL ───────────────────────────────────────────
cat(sprintf("[%s] Memory before MatrixEQTL: %.1f MB used\n",
    Sys.time(), sum(gc()[,2]) * 8 / 1024))

# ── Run MatrixEQTL ────────────────────────────────────────────────────────────
cat(sprintf("[%s] Running MatrixEQTL...\n", Sys.time()))

# MatrixEQTL accepts either a path or an open connection for its output
# arguments. A gzfile() connection compresses as it writes, so the liberal cis
# threshold never lands on disk uncompressed. When the trans pass is disabled
# (pval_trans = 0) MatrixEQTL never touches the trans output, so leave that one
# as a plain path — opening a connection there would leave a valid-but-empty
# .gz that the header-only stub below would then overwrite with plain text.
open_out <- function(path, use_gzip) {
  if (!use_gzip) return(path)
  gzfile(path, open = "wt")
}

out_cis_target   <- open_out(opt$out_cis, opt$gzip_output)
out_trans_target <- if (opt$pval_trans == 0) opt$out_trans else
                      open_out(opt$out_trans, opt$gzip_output)

me <- Matrix_eQTL_main(
  snps                  = snp_data,
  gene                  = meth_data,
  cvrt                  = covar_data,
  output_file_name      = out_trans_target,
  pvOutputThreshold     = opt$pval_trans,
  useModel              = modelLINEAR,
  errorCovariance       = numeric(),
  verbose               = TRUE,
  output_file_name.cis  = out_cis_target,
  pvOutputThreshold.cis = opt$pval_cis,
  snpspos               = snpspos_df,
  genepos               = genepos_df,
  cisDist               = opt$cis_dist,
  pvalue.hist           = FALSE,
  min.pv.by.genesnp     = FALSE,
  # Must stay FALSE. Beyond memory, noFDRsaveMemory=TRUE drops the FDR column
  # and writes results unsorted in streaming chunks. Downstream needs both: the
  # FDR column, and p-value-ascending order so the 10M-row cap that protects
  # readers from the liberal cis threshold keeps the most significant rows
  # instead of an arbitrary slice.
  noFDRsaveMemory       = FALSE
)

for (con in list(out_cis_target, out_trans_target)) {
  if (inherits(con, "connection") && isOpen(con)) close(con)
}

rm(snpspos_df, genepos_df, snp_pos, cpg_pos); gc()
cat(sprintf("[%s] Complete.\n", Sys.time()))

# When output_file_name is a connection, MatrixEQTL returns an empty me$cis —
# neqtls is NULL, and sprintf("%d", NULL) silently prints nothing. Count from
# the written file instead, streaming so a multi-GB gzip is never held in RAM.
count_rows <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con))
  n <- 0L
  repeat {
    chunk <- readLines(con, n = 1e6L, warn = FALSE)
    if (length(chunk) == 0L) break
    n <- n + length(chunk)
  }
  max(n - 1L, 0L)   # drop header
}

cat(sprintf("Cis mQTLs  (p < %.0e): %d\n", opt$pval_cis, count_rows(opt$out_cis)))
cat(sprintf("Output compression: %s\n", if (opt$gzip_output) "gzip" else "none"))

# pval_trans = 0 disables the trans pass outright: MatrixEQTL skips the
# computation (rather than merely filtering the output, which is what any
# non-zero threshold does), returns no $trans element, and never creates the
# trans output file. Snakemake declares that file as a rule output, so write a
# header-only stub to keep the rule's contract. Header matches MatrixEQTL's own
# so downstream readers need no special case.
if (is.null(me$trans)) {
  cat(sprintf("Trans mQTLs: analysis disabled (pval_trans=%g)\n", opt$pval_trans))
  stub_con <- if (opt$gzip_output) gzfile(opt$out_trans, open = "wt") else file(opt$out_trans, open = "wt")
  writeLines("SNP\tgene\tbeta\tt-stat\tp-value\tFDR", stub_con)
  close(stub_con)
  cat(sprintf("Wrote header-only trans file: %s\n", opt$out_trans))
} else {
  cat(sprintf("Trans mQTLs (p < %.0e): %d\n", opt$pval_trans, count_rows(opt$out_trans)))
}
