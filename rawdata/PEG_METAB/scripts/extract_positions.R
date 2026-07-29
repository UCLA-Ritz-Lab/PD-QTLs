#!/usr/bin/env Rscript
# extract_positions.R
# Extracts SNP positions from PLINK2 .pvar file.
# For untargeted metabolomics features, there are no genomic coordinates,
# so all features are assigned chr=0, s1=0, s2=0. This means MatrixEQTL
# treats all SNP-feature associations as trans. The cis threshold in
# config.yaml should be set equal to the trans threshold for this reason.

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--pvar_file",       type="character"),
  make_option("--feature_list",    type="character"),
  make_option("--out_snp_pos",     type="character"),
  make_option("--out_feature_pos", type="character"),
  make_option("--log",             type="character", default="extract_positions.log")
)
opt <- parse_args(OptionParser(option_list=option_list))

cat(sprintf("[%s] Extracting positions\n", Sys.time()))

# -- SNP positions from .pvar -------------------------------------------------
pvar <- fread(opt$pvar_file, skip="#CHROM", data.table=FALSE)
colnames(pvar)[1] <- "CHROM"

snp_pos <- data.frame(
  snp_id = pvar$ID,
  chr    = as.character(pvar$CHROM),
  pos    = pvar$POS
)

# Normalise chr format — strip "chr" prefix if present
snp_pos$chr <- sub("^chr", "", snp_pos$chr)

fwrite(snp_pos, opt$out_snp_pos, sep=",")
cat(sprintf("SNP positions written: %d SNPs\n", nrow(snp_pos)))

# -- Feature positions — all set to 0 for untargeted metabolomics -------------
features <- readLines(opt$feature_list)

feature_pos <- data.frame(
  geneid = features,
  chr    = "0",
  s1     = 0L,
  s2     = 0L
)

fwrite(feature_pos, opt$out_feature_pos, sep=",")
cat(sprintf("Feature positions written: %d features (all chr=0, pos=0)\n",
            nrow(feature_pos)))

cat(sprintf("[%s] extract_positions complete.\n", Sys.time()))
