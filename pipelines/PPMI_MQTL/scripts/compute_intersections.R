#!/usr/bin/env Rscript
# compute_intersections.R
# Computes the SNP and probe intersections between PEG and PPMI datasets
# for cases and controls separately. Outputs restricted position and probe
# files that are passed to run_mqtl.R to reduce memory and multiple testing.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--ppmi_snp_pos",       type="character"),
  make_option("--ppmi_probe_list",    type="character"),
  make_option("--peg_snp_pos_cases",  type="character",
              help="PEG SNP positions for cases meta-analysis"),
  make_option("--peg_snp_pos_ctrls",  type="character",
              help="PEG SNP positions for controls meta-analysis"),
  make_option("--peg_probe_cases",    type="character",
              help="PEG retained probes for cases"),
  make_option("--peg_probe_ctrls",    type="character",
              help="PEG retained probes for controls"),
  make_option("--out_dir",            type="character"),
  make_option("--log",                type="character",
              default="compute_intersections.log")
)
opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$out_dir,            recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(opt$log),       recursive=TRUE, showWarnings=FALSE)
log_con <- file(opt$log, open="wt")
sink(log_con, type="message")
sink(log_con, type="output", append=TRUE, split=TRUE)
on.exit({ sink(type="output"); sink(type="message"); close(log_con) })

log_section <- function(t) {
  cat("\n", strrep("=",60), "\n", t, "\n", strrep("=",60), "\n", sep="")
}

# ── Load PPMI positions ───────────────────────────────────────────────────────
log_section("Loading PPMI positions")

ppmi_snps <- fread(opt$ppmi_snp_pos, data.table=FALSE) %>%
  mutate(chr=as.character(chr), pos=as.integer(pos),
         chrpos=paste0(chr, ":", pos))
cat(sprintf("PPMI SNPs: %d\n", nrow(ppmi_snps)))

ppmi_probes <- readLines(opt$ppmi_probe_list)
cat(sprintf("PPMI probes: %d\n", length(ppmi_probes)))

# ── Helper: compute intersection for one meta-analysis group ──────────────────
compute_group <- function(group_name, peg_snp_file, peg_probe_file) {
  log_section(sprintf("Computing intersection: %s", group_name))

  # Load PEG positions — match by chr:pos since IDs differ
  peg_snps <- fread(peg_snp_file, data.table=FALSE) %>%
    mutate(chr=as.character(chr), pos=as.integer(pos),
           chrpos=paste0(chr, ":", pos))
  cat(sprintf("PEG SNPs: %d\n", nrow(peg_snps)))

  # SNP intersection by chr:pos
  common_chrpos <- intersect(ppmi_snps$chrpos, peg_snps$chrpos)
  cat(sprintf("Common SNP positions: %d\n", length(common_chrpos)))

  # Keep PPMI SNP entries at common positions (preserves PPMI snp_id format)
  ppmi_snps_restricted <- ppmi_snps %>%
    filter(chrpos %in% common_chrpos) %>%
    select(snp_id, chr, pos)

  cat(sprintf("PPMI SNPs retained: %d\n", nrow(ppmi_snps_restricted)))

  # Probe intersection by cg ID
  peg_probes  <- readLines(peg_probe_file)
  cat(sprintf("PEG probes: %d\n", length(peg_probes)))

  common_probes <- intersect(ppmi_probes, peg_probes)
  cat(sprintf("Common probes: %d\n", length(common_probes)))

  # Save outputs
  snp_out   <- file.path(opt$out_dir,
                          paste0(group_name, "_snp_positions.csv"))
  probe_out <- file.path(opt$out_dir,
                          paste0(group_name, "_retained_probes.txt"))

  fwrite(ppmi_snps_restricted, snp_out)
  writeLines(common_probes, probe_out)

  cat(sprintf("SNP positions saved:  %s\n", basename(snp_out)))
  cat(sprintf("Probe list saved:     %s\n", basename(probe_out)))

  invisible(list(
    n_snps   = nrow(ppmi_snps_restricted),
    n_probes = length(common_probes)
  ))
}

# ── Cases intersection ────────────────────────────────────────────────────────
result_cases <- compute_group(
  "cases",
  opt$peg_snp_pos_cases,
  opt$peg_probe_cases
)

# ── Controls intersection ─────────────────────────────────────────────────────
result_ctrls <- compute_group(
  "controls",
  opt$peg_snp_pos_ctrls,
  opt$peg_probe_ctrls
)

# ── Summary ───────────────────────────────────────────────────────────────────
log_section("Summary")
cat(sprintf("%-12s  %8s  %8s\n", "Group", "SNPs", "Probes"))
cat(strrep("-", 32), "\n")
cat(sprintf("%-12s  %8d  %8d\n", "cases",
    result_cases$n_snps, result_cases$n_probes))
cat(sprintf("%-12s  %8d  %8d\n", "controls",
    result_ctrls$n_snps, result_ctrls$n_probes))
cat(sprintf("\nEstimated MatrixEQTL matrix sizes:\n"))
cat(sprintf("  Cases:    %d SNPs x %d probes\n",
    result_cases$n_snps, result_cases$n_probes))
cat(sprintf("  Controls: %d SNPs x %d probes\n",
    result_ctrls$n_snps, result_ctrls$n_probes))
