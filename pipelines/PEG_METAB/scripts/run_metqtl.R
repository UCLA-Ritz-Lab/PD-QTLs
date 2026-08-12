#!/usr/bin/env Rscript
# run_metqtl.R
# Runs MatrixEQTL for one metQTL dataset (platform x wave x case_group).
# Modelled closely on the working PPMI run_mqtl.R for consistency.
#
# Key differences from methylation pipeline:
#   - Feature matrix is metabolomics (features x samples) not methylation
#   - Sample alignment uses pegid (metab) and GWAS_ID (geno) via id_map
#   - feature_pos has all chr=0 — all associations treated as trans
#   - SNP pre-filter still applied to reduce geno_tsv before SlicedData

suppressPackageStartupMessages({
  library(MatrixEQTL)
  library(data.table)
  library(tidyverse)
  library(optparse)
})

option_list <- list(
  make_option("--dataset",      type="character"),
  make_option("--geno_tsv",     type="character"),
  make_option("--metab_file",   type="character"),
  make_option("--covar_file",   type="character"),
  make_option("--id_map",       type="character"),
  make_option("--snp_pos",      type="character"),
  make_option("--feature_pos",  type="character"),
  make_option("--out_cis",      type="character"),
  make_option("--out_trans",    type="character"),
  make_option("--cis_dist",     type="double",  default=1e6),
  make_option("--pval_cis",     type="double",  default=1e-5),
  make_option("--pval_trans",   type="double",  default=1e-5),
  make_option("--gzip_output",  action="store_true", default=FALSE),
  make_option("--threads",      type="integer", default=1),
  make_option("--log",          type="character", default="run_metqtl.log")
)
opt <- parse_args(OptionParser(option_list=option_list))

log_con <- file(opt$log, open="wt")
sink(log_con, type="message")
sink(log_con, type="output", append=TRUE, split=TRUE)
on.exit({ sink(type="output"); sink(type="message"); close(log_con) })

cat(sprintf("[%s] Running metQTL: %s\n", Sys.time(), opt$dataset))

# ── Read genotype sample IDs from header only ─────────────────────────────────
cat(sprintf("[%s] Reading genotype header: %s\n", Sys.time(), opt$geno_tsv))
geno_header  <- fread(cmd=paste("zcat", opt$geno_tsv, "| head -1"),
                      header=FALSE, data.table=FALSE)
geno_samples <- trimws(as.character(geno_header[1, -1]))
cat(sprintf("Genotype file: %d samples in header\n", length(geno_samples)))
cat(sprintf("First 3 geno IDs: %s\n", paste(head(geno_samples, 3), collapse=", ")))
rm(geno_header); gc()

# ── Load metabolomics matrix ───────────────────────────────────────────────────
cat(sprintf("[%s] Loading metabolomics...\n", Sys.time()))
metab_df   <- fread(opt$metab_file, data.table=FALSE)
feat_ids   <- metab_df[[1]]
metab_mat  <- as.matrix(metab_df[, -1, drop=FALSE])
rownames(metab_mat) <- feat_ids
rm(metab_df); gc()

cat(sprintf("Metabolomics: %d features x %d samples\n",
    nrow(metab_mat), ncol(metab_mat)))
cat(sprintf("First 3 metab col names: %s\n",
    paste(head(colnames(metab_mat), 3), collapse=", ")))

# ── Load ID map and align samples ─────────────────────────────────────────────
# id_map columns: Sample.ID, pegid, GWAS_ID
# metabolomics columns are pegids; geno_samples are GWAS IDs
id_map <- read.csv(opt$id_map, stringsAsFactors=FALSE)
cat(sprintf("ID map: %d rows, cols: %s\n",
    nrow(id_map), paste(colnames(id_map), collapse=", ")))

metab_pegids <- trimws(colnames(metab_mat))
colnames(metab_mat) <- metab_pegids

cat(sprintf("First 3 id_map GWAS_IDs: %s\n",
    paste(head(id_map$GWAS_ID, 3), collapse=", ")))

aligned <- id_map %>%
  filter(pegid %in% metab_pegids, GWAS_ID %in% geno_samples)

cat(sprintf("Samples aligned (metab + geno): %d\n", nrow(aligned)))

if (nrow(aligned) == 0) {
  cat(sprintf("metab_pegids example: '%s'\n", metab_pegids[1]))
  cat(sprintf("geno_samples example: '%s'\n", geno_samples[1]))
  cat(sprintf("id_map pegid example: '%s'\n", id_map$pegid[1]))
  cat(sprintf("id_map GWAS_ID example: '%s'\n", id_map$GWAS_ID[1]))
  stop("No samples align between metabolomics and genotype matrices")
}

# Diagnose mismatches
pegid_in_metab <- intersect(aligned$pegid, colnames(metab_mat))
pegid_missing  <- setdiff(aligned$pegid, colnames(metab_mat))
cat(sprintf("aligned$pegid in metab_mat: %d\n", length(pegid_in_metab)))
cat(sprintf("aligned$pegid NOT in metab_mat: %d\n", length(pegid_missing)))
if (length(pegid_missing) > 0)
  aligned <- aligned %>% filter(pegid %in% colnames(metab_mat))

aligned$pegid   <- trimws(as.character(aligned$pegid))
aligned$GWAS_ID <- trimws(as.character(aligned$GWAS_ID))

metab_idx <- match(aligned$pegid,   colnames(metab_mat))
geno_idx  <- match(aligned$GWAS_ID, geno_samples)

cat(sprintf("metab_idx NAs: %d\n", sum(is.na(metab_idx))))
cat(sprintf("geno_idx NAs:  %d\n", sum(is.na(geno_idx))))

valid     <- !is.na(metab_idx) & !is.na(geno_idx)
aligned   <- aligned[valid, ]
metab_idx <- metab_idx[valid]
geno_idx  <- geno_idx[valid]

cat(sprintf("Samples after index matching: %d\n", nrow(aligned)))

metab_aligned <- metab_mat[, metab_idx, drop=FALSE]
colnames(metab_aligned) <- aligned$pegid
rm(metab_mat); gc()

# ── Load covariates ────────────────────────────────────────────────────────────
# covariates.csv written by prepare_datasets.R in vars-as-rows format:
# covariate,sample1,sample2,...
covar_df  <- read.csv(opt$covar_file, stringsAsFactors=FALSE,
                      row.names=1, check.names=FALSE)
covar_mat <- as.matrix(covar_df)
colnames(covar_mat) <- sub("^X", "", colnames(covar_mat))
rm(covar_df); gc()

cat(sprintf("Covariate matrix: %d vars x %d samples\n",
    nrow(covar_mat), ncol(covar_mat)))
cat(sprintf("First 3 covar col names: %s\n",
    paste(head(colnames(covar_mat), 3), collapse=", ")))
cat(sprintf("First 3 aligned$pegid:   %s\n",
    paste(head(aligned$pegid, 3), collapse=", ")))

# Intersect covar with aligned pegids
pegids_in_covar <- intersect(aligned$pegid, colnames(covar_mat))
pegids_missing  <- setdiff(aligned$pegid,   colnames(covar_mat))
cat(sprintf("Pegids in covar:  %d\n", length(pegids_in_covar)))
cat(sprintf("Pegids missing:   %d\n", length(pegids_missing)))

if (length(pegids_missing) > 0) {
  aligned       <- aligned %>% filter(pegid %in% pegids_in_covar)
  metab_aligned <- metab_aligned[, aligned$pegid, drop=FALSE]
}

if (nrow(aligned) == 0)
  stop("No samples remain after three-way alignment")

covar_aligned <- covar_mat[, aligned$pegid, drop=FALSE]
storage.mode(covar_aligned) <- "double"
rm(covar_mat); gc()

# Remove zero or near-zero variance covariates from the ALIGNED subset.
# A covariate that varies in the full dataset may be constant in this
# subsample (e.g. Ethnicity_Asian if all Asian samples were dropped by
# the geno alignment step). Such covariates are collinear with the
# intercept and cause MatrixEQTL to abort.
row_vars <- apply(covar_aligned, 1, var, na.rm=TRUE)
drop_covars <- names(row_vars)[is.na(row_vars) | row_vars < 1e-6]
if (length(drop_covars) > 0) {
  cat(sprintf("Dropping %d zero-variance covariates from aligned subset: %s\n",
      length(drop_covars), paste(drop_covars, collapse=", ")))
  covar_aligned <- covar_aligned[!rownames(covar_aligned) %in% drop_covars, ,
                                  drop=FALSE]
}

cat(sprintf("Final aligned samples: %d\n", ncol(covar_aligned)))
cat(sprintf("Covariates: %d vars x %d samples\n",
    nrow(covar_aligned), ncol(covar_aligned)))

# ── Load positions ─────────────────────────────────────────────────────────────
snp_pos     <- fread(opt$snp_pos,     data.table=FALSE)
feature_pos <- fread(opt$feature_pos, data.table=FALSE)

# ── Write metab and covariates to temp files ───────────────────────────────────
tmp_metab  <- tempfile(fileext=".txt")
tmp_covar  <- tempfile(fileext=".txt")
tmp_geno   <- tempfile(fileext=".tsv.gz")
tmp_snpids <- tempfile(fileext=".txt")
on.exit(unlink(c(tmp_metab, tmp_covar, tmp_geno, tmp_snpids)), add=TRUE)

write.table(metab_aligned,  tmp_metab, sep="\t", quote=FALSE)
write.table(covar_aligned,  tmp_covar, sep="\t", quote=FALSE)
rm(covar_aligned); gc()

# Pre-filter geno_tsv to SNPs in snp_pos before SlicedData loading
writeLines(snp_pos$snp_id, tmp_snpids)
cat(sprintf("[%s] Pre-filtering geno_tsv to %d SNPs...\n",
            Sys.time(), nrow(snp_pos)))

filter_cmd <- paste0(
  "{ zcat ", opt$geno_tsv, " | head -1; ",
  "zcat ", opt$geno_tsv, " | tail -n +2 | grep -Fw -f ", tmp_snpids, "; }",
  " | bgzip -c > ", tmp_geno
)
if (system(filter_cmd) != 0) stop("Failed to pre-filter geno_tsv")

n_filtered <- as.integer(system(
  paste("zcat", tmp_geno, "| tail -n +2 | wc -l"), intern=TRUE))
cat(sprintf("[%s] Filtered geno_tsv: %d SNPs retained\n", Sys.time(), n_filtered))
gc()

# ── Subset geno_tsv columns to aligned samples ────────────────────────────────
cat(sprintf("[%s] Subsetting geno_tsv to %d aligned samples...\n",
            Sys.time(), nrow(aligned)))

geno_aligned_cols <- match(aligned$GWAS_ID, geno_samples)
n_missing <- sum(is.na(geno_aligned_cols))
if (n_missing > 0) {
  keep              <- !is.na(geno_aligned_cols)
  aligned           <- aligned[keep, ]
  metab_aligned     <- metab_aligned[, aligned$pegid, drop=FALSE]
  geno_aligned_cols <- geno_aligned_cols[keep]
}

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
cat(sprintf("Metab aligned: %d samples\n", ncol(metab_aligned)))
stopifnot(n_geno_cols == ncol(metab_aligned))
rm(metab_aligned); gc()

# ── Decompress geno_aligned for SlicedData ────────────────────────────────────
tmp_geno_txt <- tempfile(fileext=".txt")
on.exit(unlink(tmp_geno_txt), add=TRUE)
if (system(paste("zcat", tmp_geno_aligned, ">", tmp_geno_txt)) != 0)
  stop("Failed to decompress geno_aligned")

# ── Build MatrixEQTL SlicedData objects ───────────────────────────────────────
snp_data    <- SlicedData$new()
metab_data  <- SlicedData$new()
covar_data  <- SlicedData$new()

snp_data$fileDelimiter      <- "\t"
snp_data$fileOmitCharacters <- "NA"
snp_data$fileSkipRows       <- 1
snp_data$fileSkipColumns    <- 1
snp_data$fileSliceSize      <- 500
snp_data$LoadFile(tmp_geno_txt)

metab_data$fileDelimiter      <- "\t"
metab_data$fileOmitCharacters <- "NA"
metab_data$fileSkipRows       <- 1
metab_data$fileSkipColumns    <- 1
metab_data$fileSliceSize      <- 2000
metab_data$LoadFile(tmp_metab)

covar_data$fileDelimiter      <- "\t"
covar_data$fileOmitCharacters <- "NA"
covar_data$fileSkipRows       <- 1
covar_data$fileSkipColumns    <- 1
covar_data$fileSliceSize      <- 1e6
covar_data$LoadFile(tmp_covar)

# ── Build position data frames ─────────────────────────────────────────────────
snpspos_df <- snp_pos %>%
  filter(snp_id %in% snp_data$GetAllRowNames()) %>%
  select(snpid=snp_id, chr, pos) %>%
  mutate(chr=as.character(chr))

genepos_df <- feature_pos %>%
  filter(geneid %in% metab_data$GetAllRowNames()) %>%
  select(geneid, chr, s1, s2) %>%
  mutate(chr=as.character(chr))

cat(sprintf("SNP positions for MatrixEQTL: %d\n", nrow(snpspos_df)))
cat(sprintf("Feature positions:            %d\n", nrow(genepos_df)))

cat(sprintf("[%s] Memory before MatrixEQTL: %.1f MB\n",
    Sys.time(), sum(gc()[,2]) * 8 / 1024))

# ── Run MatrixEQTL ─────────────────────────────────────────────────────────────
cat(sprintf("[%s] Running MatrixEQTL...\n", Sys.time()))

# MatrixEQTL writes to `output_file_name` as either a path or an open
# connection. Handing it a gzfile() connection compresses as it writes, so the
# ~8.5 GB plain trans output at pval_trans=1e-2 never lands on disk uncompressed.
open_out <- function(path, use_gzip) {
  if (!use_gzip) return(path)
  con <- gzfile(path, open = "wt")
  con
}

out_trans_target <- open_out(opt$out_trans, opt$gzip_output)
out_cis_target   <- open_out(opt$out_cis,   opt$gzip_output)

# noFDRsaveMemory must stay FALSE. It has two effects beyond memory: it drops
# the FDR column entirely, and it writes results unsorted in streaming chunks.
# Downstream consumers depend on both — the overlap report filters on FDR, and
# the row cap in config `downstream.max_trans_rows` is only meaningful because
# MatrixEQTL sorts by p-value ascending when it computes FDR, so truncating
# keeps the most significant rows rather than an arbitrary slice.
me <- Matrix_eQTL_main(
  snps                  = snp_data,
  gene                  = metab_data,
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
  noFDRsaveMemory       = FALSE
)

for (con in list(out_trans_target, out_cis_target)) {
  if (inherits(con, "connection") && isOpen(con)) close(con)
}

rm(snpspos_df, genepos_df, snp_pos, feature_pos); gc()
cat(sprintf("[%s] Complete.\n", Sys.time()))

# When output_file_name is a connection rather than a path, MatrixEQTL returns
# an empty me$cis / me$trans — neqtls is NULL and the old sprintf("%d") logging
# silently printed nothing. Count from the written file instead, streaming so a
# multi-GB gzip never has to be held in memory.
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

cat(sprintf("Cis metQTLs  (p < %.0e): %d\n", opt$pval_cis,   count_rows(opt$out_cis)))
cat(sprintf("Trans metQTLs (p < %.0e): %d\n", opt$pval_trans, count_rows(opt$out_trans)))
cat(sprintf("Output compression: %s\n", if (opt$gzip_output) "gzip" else "none"))
cat(sprintf("Trans output size: %.2f GB\n", file.size(opt$out_trans) / 1024^3))
