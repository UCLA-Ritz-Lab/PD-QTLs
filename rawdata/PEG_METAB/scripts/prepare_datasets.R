#!/usr/bin/env Rscript
# prepare_datasets.R
# Loads C18 and HILIC metabolomics RData, maps sample IDs via keyVar files,
# splits by platform x wave x case_group, filters low-variance features,
# writes per-dataset matrices, covariates, keep files, and ID maps.
#
# Covariates used:
#   Female, Age, PlasmaBlast, CD8pCD28nCD45RAn, CD4T, NK, Mono, Gran (cell types),
#   Ethnicity dummy variables (Hispanic, Asian, Black; ref=Caucasian)
#   PDstudyDiseaseNumeric is excluded — constant within case/control strata
#
# Sample ID mapping:
#   PEG1 covar: ExternalDNACode -> pegid (via keyvar$pegid)
#   PEG2 covar: PEGID           -> pegid (via keyvar$pegid)

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(optparse)
})

option_list <- list(
  make_option("--c18_rdata",    type="character"),
  make_option("--c18_keyvar",   type="character"),
  make_option("--hilic_rdata",  type="character"),
  make_option("--hilic_keyvar", type="character"),
  make_option("--covar_peg1",   type="character"),
  make_option("--covar_peg2",   type="character"),
  make_option("--gwas_linkage", type="character",
              help="CSV mapping GWAS_ID (VCF sample IDs) to Pegid"),
  make_option("--vcf",          type="character"),
  make_option("--out_dir",      type="character"),
  make_option("--min_sd",       type="double",  default=0.01),
  make_option("--max_na_frac",  type="double",  default=0.20),
  make_option("--log",          type="character", default="prepare_datasets.log")
)
opt <- parse_args(OptionParser(option_list=option_list))

log_con <- file(opt$log, open="wt")
sink(log_con, type="message")
on.exit({ sink(type="message"); close(log_con) }, add=TRUE)

cat(sprintf("[%s] Starting prepare_datasets\n", Sys.time()))

# -- Helpers ------------------------------------------------------------------
load_rdata_obj <- function(path) {
  e <- new.env()
  load(path, envir=e)
  get(ls(e)[[1]], envir=e)
}

# -- Load metabolomics data ---------------------------------------------------
c18_mat   <- as.matrix(load_rdata_obj(opt$c18_rdata))
hilic_mat <- as.matrix(load_rdata_obj(opt$hilic_rdata))
c18_key   <- load_rdata_obj(opt$c18_keyvar)
hilic_key <- load_rdata_obj(opt$hilic_keyvar)

cat(sprintf("C18:   %d samples x %d features\n", nrow(c18_mat),   ncol(c18_mat)))
cat(sprintf("HILIC: %d samples x %d features\n", nrow(hilic_mat), ncol(hilic_mat)))

# -- Load covariates ----------------------------------------------------------
covar_peg1 <- read.csv(opt$covar_peg1, stringsAsFactors=FALSE, check.names=FALSE)
covar_peg2 <- read.csv(opt$covar_peg2, stringsAsFactors=FALSE, check.names=FALSE)

# Standardise the sample ID column to "pegid" for both waves
# PEG1 uses ExternalDNACode, PEG2 uses PEGID
covar_peg1 <- covar_peg1 %>% rename(pegid = ExternalDNACode)
covar_peg2 <- covar_peg2 %>% rename(pegid = PEGID)

# Columns to use as covariates (same set for both waves)
COVAR_COLS <- c("Female", "Age",
                "PlasmaBlast", "CD8pCD28nCD45RAn", "CD4T", "NK", "Mono", "Gran")

# Dummy-code Ethnicity with Caucasian as reference level
# Produces: Ethnicity_Hispanic, Ethnicity_Asian, Ethnicity_Black
dummy_ethnicity <- function(covar_df) {
  eth <- covar_df$Ethnicity
  covar_df$Ethnicity_Hispanic <- as.integer(eth == "Hispanic")
  covar_df$Ethnicity_Asian    <- as.integer(eth == "Asian")
  covar_df$Ethnicity_Black    <- as.integer(eth == "Black")
  covar_df
}

covar_peg1 <- covar_peg1 %>%
  dummy_ethnicity() %>%
  select(pegid, all_of(COVAR_COLS),
         Ethnicity_Hispanic, Ethnicity_Asian, Ethnicity_Black)

covar_peg2 <- covar_peg2 %>%
  dummy_ethnicity() %>%
  select(pegid, all_of(COVAR_COLS),
         Ethnicity_Hispanic, Ethnicity_Asian, Ethnicity_Black)

cat(sprintf("Covar PEG1: %d samples, %d vars\n", nrow(covar_peg1), ncol(covar_peg1) - 1))
cat(sprintf("Covar PEG2: %d samples, %d vars\n", nrow(covar_peg2), ncol(covar_peg2) - 1))

# -- Get VCF sample IDs and resolve to pegids ---------------------------------
vcf_samples_raw <- system(paste("bcftools query -l", opt$vcf), intern=TRUE)
cat(sprintf("VCF samples (raw): %d, example: '%s'\n",
            length(vcf_samples_raw), vcf_samples_raw[1]))

# Load GWAS linkage file to map VCF sample IDs -> Pegid
gwas_link <- read.csv(opt$gwas_linkage, stringsAsFactors=FALSE, check.names=FALSE)
cat(sprintf("GWAS linkage: %d rows, cols: %s\n",
            nrow(gwas_link), paste(colnames(gwas_link), collapse=", ")))

# Build lookup: VCF sample ID -> Pegid
# Standard columns from PEG_WGS/gwas_linkage.csv: Pegid, GWAS_ID
if ("GWAS_ID" %in% colnames(gwas_link) && "Pegid" %in% colnames(gwas_link)) {
  gwas_to_peg <- setNames(gwas_link$Pegid, gwas_link$GWAS_ID)
} else {
  id_cols <- colnames(gwas_link)
  gwas_to_peg <- setNames(gwas_link[[id_cols[1]]], gwas_link[[id_cols[2]]])
  cat(sprintf("Using linkage columns: %s -> %s\n", id_cols[2], id_cols[1]))
}

# vcf_samples is the set of pegids with genotype data
vcf_samples <- unique(na.omit(gwas_to_peg[vcf_samples_raw]))
cat(sprintf("VCF samples mapped to pegids: %d, example: '%s'\n",
            length(vcf_samples), vcf_samples[1]))
cat(sprintf("Keyvar pegid example: '%s'\n",
    head(c18_key$pegid[!is.na(c18_key$pegid) & c18_key$pegid != "NA"], 1)))

# -- Feature filtering --------------------------------------------------------
filter_features <- function(mat, min_sd, max_na_frac, platform) {
  na_frac <- colMeans(is.na(mat))
  sd_vals  <- apply(mat, 2, sd, na.rm=TRUE)
  keep     <- na_frac <= max_na_frac & !is.na(sd_vals) & sd_vals >= min_sd
  cat(sprintf("[%s] %s: %d / %d features retained\n",
              Sys.time(), platform, sum(keep), ncol(mat)))
  mat[, keep, drop=FALSE]
}

c18_mat   <- filter_features(c18_mat,   opt$min_sd, opt$max_na_frac, "C18")
hilic_mat <- filter_features(hilic_mat, opt$min_sd, opt$max_na_frac, "HILIC")

# Prefix feature IDs with platform to avoid downstream collisions
colnames(c18_mat)   <- paste0("C18_",   colnames(c18_mat))
colnames(hilic_mat) <- paste0("HILIC_", colnames(hilic_mat))

# Save retained feature list
all_features <- c(colnames(c18_mat), colnames(hilic_mat))
dir.create(file.path(opt$out_dir, "datasets"), recursive=TRUE, showWarnings=FALSE)
writeLines(all_features,
           file.path(opt$out_dir, "datasets", "retained_features.txt"))
cat(sprintf("Total features: %d (C18: %d, HILIC: %d)\n",
            length(all_features), ncol(c18_mat), ncol(hilic_mat)))

# -- Dataset definitions ------------------------------------------------------
datasets <- list(
  c18_peg1_cases      = list(platform="c18",   wave=1, case_group="cases"),
  c18_peg1_controls   = list(platform="c18",   wave=1, case_group="controls"),
  c18_peg2_cases      = list(platform="c18",   wave=2, case_group="cases"),
  hilic_peg1_cases    = list(platform="hilic", wave=1, case_group="cases"),
  hilic_peg1_controls = list(platform="hilic", wave=1, case_group="controls"),
  hilic_peg2_cases    = list(platform="hilic", wave=2, case_group="cases")
)

# -- Process each dataset -----------------------------------------------------
for (ds_name in names(datasets)) {
  ds    <- datasets[[ds_name]]
  mat   <- if (ds$platform == "c18") c18_mat else hilic_mat
  key   <- if (ds$platform == "c18") c18_key else hilic_key
  covar <- if (ds$wave == 1) covar_peg1 else covar_peg2

  cat(sprintf("\n[%s] Processing: %s\n", Sys.time(), ds_name))

  pd_val  <- if (ds$case_group == "cases") 1 else 0

  # Strip _2 suffix from Sample.ID to match matrix rownames
  # keyvar Sample.ID: '10002AP40_2' -> matrix rowname: '10002AP40'
  key_stripped <- key %>%
    mutate(Sample.ID_stripped = sub("_[0-9]+$", "", Sample.ID))

  # Filter keyvar to this wave, case_group, and samples present in matrix + VCF
  # Also exclude string "NA" in pegid (not caught by is.na())
  key_sub <- key_stripped %>%
    filter(peg == ds$wave, PD == pd_val) %>%
    filter(!is.na(pegid), pegid != "", pegid != "NA") %>%
    filter(Sample.ID_stripped %in% rownames(mat)) %>%
    filter(pegid %in% vcf_samples)

  cat(sprintf("  After wave+PD filter: %d\n",
      sum(key_stripped$peg == ds$wave & key_stripped$PD == pd_val, na.rm=TRUE)))
  cat(sprintf("  After pegid NA filter: %d\n",
      sum(key_stripped$peg == ds$wave & key_stripped$PD == pd_val &
          !is.na(key_stripped$pegid) & key_stripped$pegid != "" &
          key_stripped$pegid != "NA", na.rm=TRUE)))
  cat(sprintf("  After Sample.ID in mat: %d\n",
      sum(key_stripped$peg == ds$wave & key_stripped$PD == pd_val &
          !is.na(key_stripped$pegid) & key_stripped$pegid != "" &
          key_stripped$pegid != "NA" &
          key_stripped$Sample.ID_stripped %in% rownames(mat), na.rm=TRUE)))
  cat(sprintf("  After pegid in VCF: %d\n", nrow(key_sub)))
  if (nrow(key_sub) == 0) {
    cat("  WARNING: No samples — check VCF sample ID format vs pegid format\n")
    cat(sprintf("  VCF sample example:   '%s'\n", vcf_samples[1]))
    cat(sprintf("  Keyvar pegid example: '%s'\n",
        head(key$pegid[!is.na(key$pegid) & key$pegid != "NA"], 1)))
    stop(sprintf("No samples for dataset %s", ds_name))
  }

  # Subset and transpose metabolomics using stripped Sample.ID
  metab_sub <- t(mat[key_sub$Sample.ID_stripped, , drop=FALSE])
  colnames(metab_sub) <- key_sub$pegid

  # Align covariates
  covar_sub <- covar %>% filter(pegid %in% key_sub$pegid)

  # Three-way intersection
  common_ids <- Reduce(intersect, list(
    key_sub$pegid, colnames(metab_sub), covar_sub$pegid
  ))
  cat(sprintf("  Final samples (three-way): %d\n", length(common_ids)))
  if (length(common_ids) == 0) stop(sprintf("No samples after alignment for %s", ds_name))

  metab_sub <- metab_sub[, common_ids, drop=FALSE]

  # Build covariate matrix: vars x samples
  covar_mat <- covar_sub %>%
    filter(pegid %in% common_ids) %>%
    arrange(match(pegid, common_ids)) %>%
    column_to_rownames("pegid") %>%
    t() %>%
    as.matrix()
  covar_mat <- covar_mat[, common_ids, drop=FALSE]
  storage.mode(covar_mat) <- "double"

  # Write outputs
  ds_dir <- file.path(opt$out_dir, "datasets", ds_name)
  dir.create(ds_dir, recursive=TRUE, showWarnings=FALSE)

  # Metabolomics: features x samples
  metab_df <- cbind(feature_id=rownames(metab_sub), as.data.frame(metab_sub))
  fwrite(metab_df, file.path(ds_dir, "metabolomics.csv.gz"), sep=",")

  # Covariates: vars x samples format matching PEG mQTL convention
  # First column is "covariate" (variable names), remaining cols are sample IDs
  # This matches the format read by run_metqtl.R with row.names=1
  covar_out <- as.data.frame(covar_mat)
  covar_out <- cbind(covariate=rownames(covar_out), covar_out)
  write.csv(covar_out, file.path(ds_dir, "covariates.csv"), row.names=FALSE)

  # ID map: pegid -> GWAS_ID (VCF sample ID from gwas_linkage)
  # gwas_to_peg was built as GWAS_ID -> pegid, so reverse it for lookup
  peg_to_gwas <- setNames(names(gwas_to_peg), gwas_to_peg)
  gwas_ids_for_common <- peg_to_gwas[common_ids]

  id_map <- data.frame(
    Sample.ID = key_sub$Sample.ID_stripped[match(common_ids, key_sub$pegid)],
    pegid     = common_ids,
    GWAS_ID   = gwas_ids_for_common
  )
  write.csv(id_map, file.path(ds_dir, "id_map.csv"), row.names=FALSE)

  # PLINK2 keep file — single column IID format, matching PEG mQTL convention
  writeLines(gwas_ids_for_common, file.path(ds_dir, "geno_keep.txt"))

  cat(sprintf("  Written: %s\n", ds_dir))
  rm(metab_sub, covar_mat, covar_out, key_sub, covar_sub); gc()
}

cat(sprintf("\n[%s] prepare_datasets complete.\n", Sys.time()))
