#!/usr/bin/env Rscript
# Same conversion logic as the convert_rds_to_csv Snakemake rule, extracted
# into a standalone script so it can run directly inside entrypoint.sh right
# after process_batch.R, instead of needing a second Batch job / sequential
# Snakemake rule for what's actually a per-batch, embarrassingly parallel step.
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--rds", type = "character"),
  make_option("--csv", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

mat <- readRDS(opt$rds)
df  <- as.data.frame(mat)
df  <- cbind(probe_id = rownames(df), df)
write.csv(df, gzfile(opt$csv), row.names = FALSE, quote = FALSE)

cat(sprintf("Converted %d probes x %d samples\n", nrow(mat), ncol(mat)))
