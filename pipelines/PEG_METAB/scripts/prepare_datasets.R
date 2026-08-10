#!/usr/bin/env Rscript
# prepare_datasets.R
# Loads C18 and HILIC metabolomics RData, maps sample IDs via keyVar files,
# splits by platform x wave x case_group, filters low-variance features,
# writes per-dataset matrices, covariates, keep files, and ID maps.
#
# Covariates used:
#   Female, Age, PlasmaBlast, CD8pCD28nCD45RAn, CD4T, NK, Mono, Gran (cell types),
#   ancestry: genetic PCs from --ancestry_pcs when supplied, otherwise the
#     legacy self-reported Ethnicity dummies (Hispanic, Asian, Black;
#     ref=Caucasian). The two are alternatives, never both — self-report and
#     the leading PCs encode the same axis and would be near-collinear.
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
  make_option("--ancestry_pcs", type="character", default=NULL),
  make_option("--n_pcs",        type="integer",   default=5L),
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

# Ethnicity is carried through unencoded for now; the ancestry section below
# either drops it in favour of genetic PCs or dummy-codes it as the fallback.
covar_peg1 <- covar_peg1 %>% select(pegid, all_of(COVAR_COLS), Ethnicity)
covar_peg2 <- covar_peg2 %>% select(pegid, all_of(COVAR_COLS), Ethnicity)

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

# -- Ancestry adjustment ------------------------------------------------------
# Genetic PCs are keyed by GWAS_ID (CRG_*/CRG2_*), the VCF's sample namespace,
# while every covariate/metabolomics table here is keyed by pegid. They must be
# joined through the linkage — intersecting the two namespaces directly gives
# the empty set. Uses the DUP=0 rows only, so one pegid maps to one genotype.
use_pcs <- !is.null(opt$ancestry_pcs) && nzchar(opt$ancestry_pcs)

if (use_pcs) {
  if (!all(c("Pegid", "GWAS_ID", "DUP") %in% colnames(gwas_link))) {
    stop(sprintf(paste0(
      "Ancestry PC join needs Pegid, GWAS_ID and DUP columns in %s.\n",
      "  Found: %s"),
      opt$gwas_linkage, paste(colnames(gwas_link), collapse=", ")))
  }

  linkage_dedup <- gwas_link %>%
    filter(DUP == 0, !is.na(Pegid), Pegid != "", !is.na(GWAS_ID), GWAS_ID != "") %>%
    select(Pegid, GWAS_ID)
  cat(sprintf("Linkage rows (DUP=0): %d\n", nrow(linkage_dedup)))

  anc_raw <- read.csv(opt$ancestry_pcs, stringsAsFactors=FALSE)
  avail   <- grep("^PC[0-9]+$", colnames(anc_raw), value=TRUE)
  avail   <- avail[order(as.integer(sub("^PC", "", avail)))]
  if (opt$n_pcs > length(avail)) {
    stop(sprintf("Requested %d ancestry PCs but file has only %d (%s)",
                 opt$n_pcs, length(avail), opt$ancestry_pcs))
  }
  pc_cols <- head(avail, opt$n_pcs)

  anc_pc <- anc_raw %>%
    select(GWAS_ID = IID, all_of(pc_cols)) %>%
    inner_join(linkage_dedup, by="GWAS_ID") %>%
    select(pegid = Pegid, all_of(pc_cols))

  cat(sprintf("Ancestry PC file: %d samples, %d PCs available, using %s\n",
              nrow(anc_raw), length(avail), paste(pc_cols, collapse=", ")))
  cat(sprintf("Resolved to pegid via linkage: %d\n", nrow(anc_pc)))
  if (nrow(anc_pc) == 0) {
    stop(sprintf(paste0(
      "No ancestry PC sample resolved to a pegid.\n",
      "  Example PC IID:          %s\n",
      "  Example linkage GWAS_ID: %s"),
      paste(head(anc_raw$IID, 3),           collapse=", "),
      paste(head(linkage_dedup$GWAS_ID, 3), collapse=", ")))
  }
  if (anyDuplicated(anc_pc$pegid)) {
    stop("Ancestry PCs resolved to duplicate pegids — check DUP filtering in linkage")
  }

  # Ethnicity dummies are dropped here, not kept alongside: the leading PCs and
  # the self-report encode the same axis and would be near-collinear.
  #
  # inner_join, not left_join: a sample with no PCs would otherwise carry NA
  # covariates into MatrixEQTL. Dropping it from the covariate table lets the
  # existing three-way intersection exclude it, which is where every other
  # sample loss in this script is already accounted for.
  #
  # PCs are standardised because projected scores run from ~0.15 (PC1) down to
  # ~1e-4 (PC5); left raw, the small-scale PCs yield correspondingly huge betas.
  # This is an affine rescaling, so MatrixEQTL p-values are unchanged. Done per
  # wave — also affine, so per-dataset would give identical results.
  attach_pcs <- function(covar_df, wave_label) {
    n_before <- nrow(covar_df)
    out <- covar_df %>%
      select(-Ethnicity) %>%
      inner_join(anc_pc, by="pegid")
    cat(sprintf("%s: %d / %d samples have ancestry PCs\n",
                wave_label, nrow(out), n_before))
    if (nrow(out) == 0) {
      stop(sprintf("No %s sample matched the ancestry PC file", wave_label))
    }
    for (pc in pc_cols) {
      mu  <- mean(out[[pc]], na.rm=TRUE)
      sdv <- sd(out[[pc]],   na.rm=TRUE)
      if (is.na(sdv) || sdv == 0) {
        cat(sprintf("  WARNING: %s has zero variance in %s — left raw\n",
                    pc, wave_label))
        next
      }
      out[[pc]] <- (out[[pc]] - mu) / sdv
    }
    out
  }

  covar_peg1 <- attach_pcs(covar_peg1, "PEG1")
  covar_peg2 <- attach_pcs(covar_peg2, "PEG2")
  cat(sprintf("Ancestry PCs standardised per wave: %s\n",
              paste(pc_cols, collapse=", ")))
} else {
  cat("No ancestry PC file supplied — falling back to self-reported Ethnicity dummies\n")
  covar_peg1 <- dummy_ethnicity(covar_peg1) %>% select(-Ethnicity)
  covar_peg2 <- dummy_ethnicity(covar_peg2) %>% select(-Ethnicity)
}

cat(sprintf("Covar PEG1: %d samples, %d vars (%s)\n",
            nrow(covar_peg1), ncol(covar_peg1) - 1,
            paste(setdiff(colnames(covar_peg1), "pegid"), collapse=", ")))
cat(sprintf("Covar PEG2: %d samples, %d vars (%s)\n",
            nrow(covar_peg2), ncol(covar_peg2) - 1,
            paste(setdiff(colnames(covar_peg2), "pegid"), collapse=", ")))

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
