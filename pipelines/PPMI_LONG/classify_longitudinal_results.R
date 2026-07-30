library(dplyr)
library(readr)

# ============================================================
# Step 1-2: FDR-correct longitudinal mixed-model results and
# classify each (SNP, CpG) pair against the meta-analysis hit list.
#
# Run separately per group (cases/controls) since they are
# distinct populations with their own multiple-testing burden.
# ============================================================

classify_group <- function(group, meta_file) {

  cat("\n========== Classifying:", group, "==========\n")

  results <- read_csv(sprintf("../../results/PPMI_LONG/longitudinal_mqtl_%s_results.csv", group),
                       col_types = cols(pair_id = col_character(), snp_id = col_character(),
                                        cpg_id = col_character(), status = col_character(),
                                        warnings = col_character(), .default = col_double()),
                       guess_max = Inf)

  cat("Total pairs tested:", nrow(results), "\n")
  cat("Status breakdown:\n")
  print(table(results$status))

  # FDR correction only makes sense over pairs that actually produced a model —
  # skipped_low_n / model_failed rows have no p-value to correct and are
  # excluded from the p.adjust() denominator (not just from the results table).
  tested <- results %>% filter(status == "ok")
  cat("Pairs entering FDR correction:", nrow(tested), "\n")

  tested <- tested %>%
    mutate(
      main_fdr        = p.adjust(main_p, method = "BH"),
      interaction_fdr = p.adjust(interaction_p, method = "BH")
    )

  # Bring in the original meta-analysis effect for sign concordance
  meta_hits <- read_csv(meta_file, col_types = cols(.default = "c"), guess_max = Inf) %>%
    mutate(meta_beta = as.numeric(meta_beta), meta_pval = as.numeric(meta_pval),
           meta_FDR_BH = as.numeric(FDR_BH)) %>%
    select(pair_id, meta_beta, meta_pval, meta_FDR_BH)

  classified <- tested %>%
    inner_join(meta_hits, by = "pair_id") %>%
    mutate(
      same_sign_as_meta = sign(main_beta) == sign(meta_beta),
      replicated = main_fdr < 0.05 & same_sign_as_meta,
      dynamic    = interaction_fdr < 0.05,
      category = case_when(
        replicated & dynamic  ~ "replicated_and_dynamic",
        replicated & !dynamic ~ "replicated_static",
        !replicated & dynamic ~ "dynamic_only",       # interaction signif despite main effect not confirming
        TRUE                  ~ "not_replicated"
      )
    )

  cat("\nMain-effect FDR < 0.05:", sum(classified$main_fdr < 0.05), "\n")
  cat("Sign-concordant with meta:", sum(classified$same_sign_as_meta), "/", nrow(classified), "\n")
  cat("Interaction FDR < 0.05 (dynamic):", sum(classified$dynamic), "\n")
  cat("\nCategory breakdown:\n")
  print(table(classified$category))

  out_file <- sprintf("../../results/PPMI_LONG/longitudinal_mqtl_%s_classified.csv", group)
  write_csv(classified, out_file)
  cat("Written:", out_file, "\n")

  classified
}

cases_classified    <- classify_group("cases",    "../../results/META_MQTL/meta/cases_sig_fdr05.csv")
controls_classified <- classify_group("controls", "../../results/META_MQTL/meta/controls_sig_fdr05.csv")

cat("\n========== DONE ==========\n")
