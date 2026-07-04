#!/usr/bin/env Rscript
# extract_positions.R
# Extracts SNP positions from PLINK2 .pvar and CpG positions
# from the EPIC v1 annotation for MatrixEQTL cis/trans classification.
# Reference genome: GRCh37/hg19 (PPMI WGS pipeline reference)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(minfi)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
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

dir.create(dirname(opt$log), recursive=TRUE, showWarnings=FALSE)
log_con <- file(opt$log, open="wt")
sink(log_con, type="message")
sink(log_con, type="output", append=TRUE, split=TRUE)
on.exit({ sink(type="output"); sink(type="message"); close(log_con) })

# ── SNP positions from PLINK2 .pvar ──────────────────────────────────────────
cat("Loading SNP positions from .pvar...\n")
# Count ## meta-header lines to skip (VCF-style headers in pvar)
# fread skip= can take a string to search for, or an integer line count
con        <- file(opt$pvar_file, "r")
skip_lines <- 0L
repeat {
  ln <- readLines(con, n=1)
  if (length(ln) == 0 || !startsWith(ln, "##")) break
  skip_lines <- skip_lines + 1L
}
close(con)
cat(sprintf("Skipping %d ## header lines\n", skip_lines))

pvar_raw <- fread(
  opt$pvar_file,
  data.table = FALSE,
  skip       = skip_lines
)
cat(sprintf("pvar columns: %s\n", paste(colnames(pvar_raw), collapse=", ")))

# --set-missing-var-ids ensures ID column exists in format chr:pos:ref:alt
# matching the .raw export format from export_genotypes step
if (!"ID" %in% colnames(pvar_raw)) {
  stop(paste("No ID column in pvar — check --set-missing-var-ids was used.",
             "Columns found:", paste(colnames(pvar_raw), collapse=", ")))
}

pvar <- pvar_raw %>%
  select(snp_id=ID, chr=`#CHROM`, pos=POS) %>%
  mutate(chr=as.character(chr))

cat(sprintf("SNPs: %d\n", nrow(pvar)))
cat(sprintf("Example SNP IDs: %s\n", paste(head(pvar$snp_id, 3), collapse=", ")))
fwrite(pvar, opt$out_snp_pos)
cat(sprintf("SNP positions saved: %s\n", opt$out_snp_pos))

# ── CpG positions from EPIC v1 annotation (hg19) ──────────────────────────────
cat("Loading CpG positions from EPIC v1 annotation (hg19)...\n")

anno <- as.data.frame(
  getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
)

probe_list <- readLines(opt$probe_list)

cpg_pos <- anno %>%
  filter(Name %in% probe_list) %>%
  select(cpg_id=Name, chr=chr, pos=pos) %>%
  mutate(chr=str_remove(chr, "^chr"))   # strip chr prefix — PPMI uses bare 1,2,3

n_annotated <- nrow(cpg_pos)
n_missing   <- length(probe_list) - n_annotated
cat(sprintf("CpGs annotated: %d / %d\n", n_annotated, length(probe_list)))
if (n_missing > 0) {
  cat(sprintf("WARNING: %d probes not found in EPIC v1 annotation\n", n_missing))
}

fwrite(cpg_pos, opt$out_cpg_pos)
cat(sprintf("CpG positions saved: %s\n", opt$out_cpg_pos))
