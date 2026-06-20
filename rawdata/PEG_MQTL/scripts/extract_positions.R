#!/usr/bin/env Rscript
# extract_positions.R
# Extracts SNP positions from PLINK2 .pvar and CpG positions
# from the 450K annotation for MatrixEQTL cis/trans classification.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(optparse)
})

option_list <- list(
  make_option("--pvar_file",   type="character"),
  make_option("--probe_list",  type="character"),
  make_option("--out_snp_pos", type="character"),
  make_option("--out_cpg_pos", type="character"),
  make_option("--log",         type="character", default="extract_positions.log")
)
opt <- parse_args(OptionParser(option_list=option_list))

log_con <- file(opt$log, open="wt")
sink(log_con, type="message")
on.exit({ sink(type="message"); close(log_con) })

# ── SNP positions from PLINK2 .pvar ──────────────────────────────────────────
cat("Loading SNP positions from .pvar...\n")
pvar <- fread(opt$pvar_file, data.table=FALSE) %>%
  select(snp_id=ID, chr=`#CHROM`, pos=POS) %>%
  mutate(chr=as.character(chr))

cat(sprintf("SNPs: %d\n", nrow(pvar)))
fwrite(pvar, opt$out_snp_pos)
cat(sprintf("SNP positions saved: %s\n", opt$out_snp_pos))

# ── CpG positions from 450K annotation ───────────────────────────────────────
cat("Loading CpG positions from 450K annotation...\n")
anno <- as.data.frame(
  getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
)

probe_list <- readLines(opt$probe_list)

cpg_pos <- anno %>%
  filter(Name %in% probe_list) %>%
  select(cpg_id=Name, chr=chr, pos=pos) %>%
  mutate(chr=str_remove(chr, "^chr"))   # strip chr prefix for consistency

cat(sprintf("CpGs annotated: %d / %d requested\n",
    nrow(cpg_pos), length(probe_list)))

fwrite(cpg_pos, opt$out_cpg_pos)
cat(sprintf("CpG positions saved: %s\n", opt$out_cpg_pos))
