library(dplyr)
library(readr)

# ============================================================
# Step 1: Combine c18 + hilic PEG1 trans-metQTL results per group,
#         BH-correct once across the combined set, filter to FDR-sig.
# Step 2: Overlap with classified longitudinal mQTL results (by SNP).
#
# NOTE: metabolites carry no real genomic coordinate (chr0:pos0
# placeholder), so "trans" here just means genome-wide SNP testing,
# not a true cis/trans split. PEG2 intentionally excluded (kept as
# a sensitivity set later), mirroring the mQTL meta-analysis design.
# ============================================================

combine_and_fdr <- function(group_label, files, fdr_threshold = 0.05) {
  cat("\n---", group_label, "---\n")

  combined <- bind_rows(lapply(files, function(f) {
    read_tsv(f, col_types = cols(SNP = col_character(), gene = col_character(),
                                  beta = col_double(), `t-stat` = col_double(),
                                  `p-value` = col_double(), FDR = col_double()),
              guess_max = Inf) %>%
      mutate(source_file = basename(dirname(f)))
  }))

  cat("Total tests combined:", nrow(combined), "\n")
  cat("From files:\n"); print(table(combined$source_file))

  # NOTE: trans_metqtls.txt files were generated with a liberal pvOutputThreshold
  # (1e-2) AND manually truncated to the first 10M p-value-sorted rows to avoid
  # OOM on load. This means the file is NOT the full test space, so a BH
  # correction computed fresh on this file would be invalid (it would rank
  # already-preselected low p-values against each other rather than against
  # the true, much larger, denominator of all tests performed).
  #
  # MatrixEQTL's own FDR column, by contrast, IS valid — it was computed
  # internally using the true total test count before either truncation step.
  # Confirmed safe: FDR at the last (10-millionth) row of every source file is
  # ~0.75-0.85, well above 0.05, meaning truncation only cut into the
  # non-significant tail and nothing FDR<0.05-significant was lost.
  #
  # Consequence: c18 and hilic are treated as separate multiple-testing
  # families (each file's FDR reflects only that platform's own test space),
  # not one combined correction — not recoverable from truncated data anyway.
  n_sig <- sum(combined$FDR < fdr_threshold)
  cat("FDR <", fdr_threshold, "(per-file MatrixEQTL FDR):", n_sig, "/", nrow(combined), "\n")

  sig <- combined %>% filter(FDR < fdr_threshold)

  out_all <- sprintf("peg1_%s_metqtl_combined_fdr.csv", group_label)
  out_sig <- sprintf("peg1_%s_metqtl_sig_fdr05.csv", group_label)
  write_csv(combined, out_all)
  write_csv(sig, out_sig)
  cat("Written:", out_all, "and", out_sig, "\n")

  sig
}

metqtl_cases_sig <- combine_and_fdr(
  "cases",
  c("../PEG_METAB/results/metqtl/c18_peg1_cases/trans_metqtls.txt",
    "../PEG_METAB/results/metqtl/hilic_peg1_cases/trans_metqtls.txt")
)

metqtl_controls_sig <- combine_and_fdr(
  "controls",
  c("../PEG_METAB/results/metqtl/c18_peg1_controls/trans_metqtls.txt",
    "../PEG_METAB/results/metqtl/hilic_peg1_controls/trans_metqtls.txt")
)

# ============================================================
# Overlap with classified longitudinal mQTL results, by SNP,
# cases-vs-cases and controls-vs-controls kept separate.
# ============================================================

overlap_group <- function(group_label, mqtl_classified_file, metqtl_sig) {
  cat("\n=== Overlap:", group_label, "===\n")

  mqtl <- read_csv(mqtl_classified_file,
                    col_types = cols(pair_id = col_character(), snp_id = col_character(),
                                      cpg_id = col_character(), status = col_character(),
                                      warnings = col_character(), category = col_character(),
                                      same_sign_as_meta = col_logical(), replicated = col_logical(),
                                      dynamic = col_logical(), .default = col_double()),
                    guess_max = Inf)

  cat("mQTL classified pairs:", nrow(mqtl), "\n")
  cat("metQTL FDR-sig SNP-metabolite pairs:", nrow(metqtl_sig), "\n")

  overlap <- mqtl %>%
    inner_join(metqtl_sig, by = c("snp_id" = "SNP")) %>%
    rename(metabolite_feature = gene, metqtl_beta = beta, metqtl_pval = `p-value`, metqtl_fdr = FDR)

  cat("Unique SNPs shared between mQTL hits and metQTL hits:",
      n_distinct(overlap$snp_id), "\n")
  cat("Total overlap rows (SNP x CpG x metabolite combinations):", nrow(overlap), "\n")
  cat("Category breakdown among overlapping SNPs:\n")
  print(table(overlap %>% distinct(snp_id, category) %>% pull(category)))

  out_file <- sprintf("mqtl_metqtl_overlap_%s.csv", group_label)
  write_csv(overlap, out_file)
  cat("Written:", out_file, "\n")

  overlap
}

overlap_cases <- overlap_group("cases",
  "../PPMI_LONG/longitudinal_mqtl_cases_classified.csv", metqtl_cases_sig)

overlap_controls <- overlap_group("controls",
  "../PPMI_LONG/longitudinal_mqtl_controls_classified.csv", metqtl_controls_sig)

cat("\n========== DONE ==========\n")
