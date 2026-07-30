#!/usr/bin/env Rscript
# prepare_datasets.R
# Prepares analysis-ready methylation and covariate matrices for all
# PEG datasets. Called once by the prepare_datasets Snakemake rule.
# Outputs per-dataset files to results/datasets/{dataset}/

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--meth_peg1",        type="character"),
  make_option("--meth_peg2",        type="character"),
  make_option("--sentrix_map_peg1", type="character"),
  make_option("--sentrix_map_peg2", type="character"),
  make_option("--covar_peg1",       type="character"),
  make_option("--covar_peg2",       type="character"),
  make_option("--gwas_linkage",     type="character"),
  make_option("--vcf",              type="character"),
  make_option("--out_dir",          type="character",
             help="Root results directory — datasets/ and logs/ created here"),
  make_option("--min_sd",           type="double",  default=0.02),
  make_option("--max_na_frac",      type="double",  default=0.00),
  make_option("--log",              type="character", default="prepare_datasets.log")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Create log directory if needed and start logging immediately
dir.create(dirname(opt$log), recursive=TRUE, showWarnings=FALSE)
log_con <- file(opt$log, open="wt")
sink(log_con, type="message")
# Also echo stdout to log
sink(log_con, type="output", append=TRUE, split=TRUE)
on.exit({
  sink(type="output")
  sink(type="message")
  close(log_con)
})

log_section <- function(t) {
  cat("\n", strrep("=",60), "\n", t, "\n", strrep("=",60), "\n", sep="")
}

# ── Load linkage and VCF samples ──────────────────────────────────────────────
log_section("Step 1: ID linkage")

linkage <- read.csv(opt$gwas_linkage, stringsAsFactors=FALSE) %>%
  filter(DUP==0, Pegid!="", GWAS_ID!="") %>%
  select(Pegid, GWAS_ID)
cat(sprintf("Linkage rows (DUP=0): %d\n", nrow(linkage)))

vcf_samples <- system(sprintf("bcftools query -l %s", opt$vcf), intern=TRUE)
cat(sprintf("VCF samples: %d\n", length(vcf_samples)))
cat(sprintf("GWAS_IDs in VCF: %d / %d\n",
    sum(linkage$GWAS_ID %in% vcf_samples), nrow(linkage)))

# ── Load PEG1 ─────────────────────────────────────────────────────────────────
log_section("Step 2: PEG1")

load(opt$meth_peg1)
peg1_mat <- t(datMethPEG1t)
peg1_mat <- peg1_mat[grepl("^cg", rownames(peg1_mat)), ]
colnames(peg1_mat) <- sub("^X", "", colnames(peg1_mat))
cat(sprintf("PEG1: %d probes x %d samples\n", nrow(peg1_mat), ncol(peg1_mat)))

sentrix1 <- read.csv(opt$sentrix_map_peg1, stringsAsFactors=FALSE) %>%
  mutate(SampleID = sub("^X","",SampleID)) %>%
  rename(sentrix_id=SampleID, Pegid=ExternalDNACode)

covar1 <- read.csv(opt$covar_peg1, stringsAsFactors=FALSE) %>%
  rename(Pegid=ExternalDNACode)

peg1_meta <- data.frame(sentrix_id=colnames(peg1_mat), stringsAsFactors=FALSE) %>%
  left_join(sentrix1, by="sentrix_id") %>%
  left_join(covar1,   by="Pegid") %>%
  mutate(
    wave         = "PEG1",
    case_control = case_when(
      PDstudyDiseaseNumeric==1 ~ "case",
      PDstudyDiseaseNumeric==0 ~ "control",
      TRUE ~ NA_character_
    )
  )

cat(sprintf("PEG1 cases:    %d\n", sum(peg1_meta$case_control=="case",    na.rm=TRUE)))
cat(sprintf("PEG1 controls: %d\n", sum(peg1_meta$case_control=="control", na.rm=TRUE)))
cat(sprintf("PEG1 missing:  %d\n", sum(is.na(peg1_meta$case_control))))

valid_peg1 <- peg1_meta %>% filter(!is.na(Pegid), !is.na(PDstudyDiseaseNumeric))
peg1_mat_pid <- peg1_mat[, valid_peg1$sentrix_id, drop=FALSE]
colnames(peg1_mat_pid) <- valid_peg1$Pegid

# ── Load PEG2 ─────────────────────────────────────────────────────────────────
log_section("Step 3: PEG2")

load(opt$meth_peg2)
peg2_mat <- t(datMethPEG2t)
peg2_mat <- peg2_mat[grepl("^cg", rownames(peg2_mat)), ]
if (storage.mode(peg2_mat)=="character") {
  peg2_mat <- matrix(as.numeric(peg2_mat), nrow=nrow(peg2_mat),
                      ncol=ncol(peg2_mat), dimnames=dimnames(peg2_mat))
}
colnames(peg2_mat) <- sub("^X","",colnames(peg2_mat))
cat(sprintf("PEG2: %d probes x %d samples\n", nrow(peg2_mat), ncol(peg2_mat)))

sentrix2 <- read.csv(opt$sentrix_map_peg2, stringsAsFactors=FALSE) %>%
  mutate(SampleID=sub("^X","",SampleID)) %>%
  rename(sentrix_id=SampleID, Pegid=PEGID)

covar2 <- read.csv(opt$covar_peg2, stringsAsFactors=FALSE) %>%
  rename(Pegid=PEGID)

peg2_meta <- data.frame(sentrix_id=colnames(peg2_mat), stringsAsFactors=FALSE) %>%
  left_join(sentrix2, by="sentrix_id") %>%
  left_join(covar2,   by="Pegid") %>%
  mutate(wave="PEG2", case_control="case")

unmatched2 <- peg2_meta %>% filter(is.na(Pegid))
if (nrow(unmatched2)>0) {
  cat(sprintf("PEG2 unmatched Sentrix IDs (%d):\n", nrow(unmatched2)))
  print(unmatched2$sentrix_id)
}

valid_peg2 <- peg2_meta %>% filter(!is.na(Pegid))
peg2_mat_pid <- peg2_mat[, valid_peg2$sentrix_id, drop=FALSE]
colnames(peg2_mat_pid) <- valid_peg2$Pegid

# ── Probe filtering ───────────────────────────────────────────────────────────
log_section("Step 4: Probe filtering")

common_probes <- intersect(rownames(peg1_mat_pid), rownames(peg2_mat_pid))
cat(sprintf("Common probes: %d\n", length(common_probes)))

probe_sds <- apply(peg1_mat_pid[common_probes,], 1, sd, na.rm=TRUE)
probe_nas <- rowMeans(is.na(peg1_mat_pid[common_probes,]))

keep_probes <- names(probe_sds)[
  probe_sds >= opt$min_sd & probe_nas <= opt$max_na_frac
]
cat(sprintf("Removed (SD < %.2f):      %d\n", opt$min_sd, sum(probe_sds < opt$min_sd)))
cat(sprintf("Removed (NA > %.2f):      %d\n", opt$max_na_frac,
    sum(probe_nas > opt$max_na_frac & probe_sds >= opt$min_sd)))
cat(sprintf("Probes retained:           %d\n", length(keep_probes)))

peg1_filt <- peg1_mat_pid[keep_probes,]
peg2_filt <- peg2_mat_pid[keep_probes,]

# Save retained probe list
writeLines(keep_probes, file.path(opt$out_dir, "datasets", "retained_probes.txt"))

# ── Covariate columns ─────────────────────────────────────────────────────────
covar_cols <- c("Female","Age","PlasmaBlast","CD8pCD28nCD45RAn",
                "CD4T","NK","Mono","Gran")

make_covar_mat <- function(meta_df, cols) {
  cols     <- intersect(cols, colnames(meta_df))
  sub_df   <- meta_df %>%
    select(Pegid, all_of(cols)) %>%
    filter(!is.na(Pegid))
  # Build covariates x samples matrix (MatrixEQTL convention)
  mat           <- t(as.matrix(sub_df[, cols, drop=FALSE]))
  colnames(mat) <- sub_df$Pegid
  rownames(mat) <- cols
  # Report any NA values in covariates
  na_counts <- rowSums(is.na(mat))
  if (any(na_counts > 0)) {
    cat("  WARNING: NA values in covariates:
")
    print(na_counts[na_counts > 0])
  }
  mat
}

peg1_covar <- make_covar_mat(valid_peg1, covar_cols)
peg2_covar <- make_covar_mat(valid_peg2, covar_cols)

# ── Save per-dataset outputs ──────────────────────────────────────────────────
log_section("Step 5: Save datasets")

save_dataset <- function(beta_mat, covar_mat, meta_df,
                          dataset_name, gwas_lnk, vcf_samps) {
  out <- file.path(opt$out_dir, "datasets", dataset_name)
  dir.create(out, recursive=TRUE, showWarnings=FALSE)

  pegids   <- colnames(beta_mat)
  gwas_map <- gwas_lnk %>%
    filter(Pegid %in% pegids, GWAS_ID %in% vcf_samps)

  cat(sprintf("\n[%s]\n", dataset_name))
  cat(sprintf("  Methylation samples:  %d\n", length(pegids)))
  cat(sprintf("  Matched to VCF:       %d\n", nrow(gwas_map)))

  if (nrow(gwas_map)==0) {
    cat("  WARNING: no VCF matches — skipping\n")
    return(invisible(NULL))
  }

  matched <- gwas_map$Pegid
  b_sub   <- beta_mat[, matched, drop=FALSE]
  c_sub   <- covar_mat[, matched, drop=FALSE]
  m_sub   <- meta_df %>% filter(Pegid %in% matched)

  # Methylation (probe_id + Pegid columns)
  mf <- file.path(out, "methylation.csv.gz")
  fwrite(cbind(probe_id=rownames(b_sub), as.data.frame(b_sub)), mf)

  # Covariates (covariate + Pegid columns)
  cf <- file.path(out, "covariates.csv")
  write.csv(cbind(covariate=rownames(c_sub), as.data.frame(c_sub)),
            cf, row.names=FALSE, quote=FALSE)

  # PLINK2 --keep file (GWAS_ID, single column)
  # The PEG VCF has no FID column in its psam (#IID SEX only).
  # PLINK2 accepts a single-column --keep file when there is no FID.
  gf <- file.path(out, "geno_keep.txt")
  writeLines(gwas_map$GWAS_ID, gf)

  # Pegid → GWAS_ID map
  write.csv(gwas_map, file.path(out, "id_map.csv"),
            row.names=FALSE, quote=FALSE)

  # Sample metadata
  write.csv(m_sub %>% left_join(gwas_map, by="Pegid"),
            file.path(out, "sample_metadata.csv"),
            row.names=FALSE, quote=FALSE)

  cat(sprintf("  Saved to %s\n", out))
  invisible(list(n_samples=nrow(gwas_map), n_probes=nrow(b_sub)))
}

peg1_case_ids <- valid_peg1 %>% filter(case_control=="case") %>% pull(Pegid)
peg1_ctrl_ids <- valid_peg1 %>% filter(case_control=="control") %>% pull(Pegid)

save_dataset(peg1_filt[,peg1_case_ids], peg1_covar[,peg1_case_ids],
             valid_peg1 %>% filter(case_control=="case"),
             "peg1_cases", linkage, vcf_samples)

save_dataset(peg1_filt[,peg1_ctrl_ids], peg1_covar[,peg1_ctrl_ids],
             valid_peg1 %>% filter(case_control=="control"),
             "peg1_controls", linkage, vcf_samples)

save_dataset(peg2_filt, peg2_covar, valid_peg2,
             "peg2_cases", linkage, vcf_samples)

log_section("Done")
cat("Output directory:", file.path(opt$out_dir, "datasets"), "\n")
