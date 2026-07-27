library(dplyr)
library(readr)
library(stringr)
library(purrr)

# ============================================================
# Standalone re-run of Step 3 only. The expensive annotation step
# already completed (5/6 db x platform combos have full Stage1-5
# output on disk; c18/LipidMaps stopped at Stage3, as expected given
# it only had 2 initial DB matches). This just re-derives the
# consolidated output files from what's already there.
# ============================================================

overlap_cases    <- read_csv("mqtl_metqtl_overlap_cases_annotated.csv",
                              col_types = cols(pair_id = col_character(), snp_id = col_character(),
                                                cpg_id = col_character(), status = col_character(),
                                                warnings = col_character(), category = col_character(),
                                                metabolite_feature = col_character(),
                                                nearest_gene = col_character(), gene_relationship = col_character(),
                                                same_sign_as_meta = col_logical(), replicated = col_logical(),
                                                dynamic = col_logical(), .default = col_double()),
                              guess_max = Inf)
overlap_controls <- read_csv("mqtl_metqtl_overlap_controls_annotated.csv",
                              col_types = cols(pair_id = col_character(), snp_id = col_character(),
                                                cpg_id = col_character(), status = col_character(),
                                                warnings = col_character(), category = col_character(),
                                                metabolite_feature = col_character(),
                                                nearest_gene = col_character(), gene_relationship = col_character(),
                                                same_sign_as_meta = col_logical(), replicated = col_logical(),
                                                dynamic = col_logical(), .default = col_double()),
                              guess_max = Inf)

parse_feature <- function(feature) {
  platform <- str_extract(feature, "^[A-Za-z0-9]+")
  nums <- str_extract(feature, "(?<=mz_rt_).*") %>% str_split("_")
  tibble(
    metabolite_feature = feature,
    platform = platform,
    mz = as.numeric(map_chr(nums, 1)),
    rt = as.numeric(map_chr(nums, 2))
  )
}

overlap_features <- bind_rows(overlap_cases, overlap_controls) %>%
  distinct(metabolite_feature) %>%
  filter(!is.na(metabolite_feature)) %>%
  pull(metabolite_feature) %>%
  map_dfr(parse_feature)

cat("Overlap features to match:", nrow(overlap_features), "\n")

MASS_TOL_PPM <- 10
RT_TOL_SEC <- 30

ppm <- function(theo_mz, obs_mz) ((obs_mz - theo_mz) / theo_mz) * 1e6

match_annotator <- function(id, mzr, rt, stage5_tbl) {
  tbl <- stage5_tbl %>%
    filter(abs(ppm(round(mz, 3), round(mzr, 3))) <= MASS_TOL_PPM) %>%
    filter(abs(time - rt) <= RT_TOL_SEC) %>%
    rename(match_chemical = name, rt = time) %>%
    select(mz, rt, match_chemical, chemical_id, confidence, score)

  if (nrow(tbl) == 0) {
    tibble(id = id, mz = mzr, rt = rt, match_chemical = "",
           chemical_id = "", confidence = NA_real_, score = NA_real_)
  } else {
    tbl %>% mutate(id = id) %>% relocate(id)
  }
}

stage5_files <- list.files("annotation_overlap", pattern = "^Stage5\\.csv$",
                            full.names = TRUE, recursive = TRUE)

cat("Stage5.csv files found:", length(stage5_files), "\n")
print(stage5_files)

if (length(stage5_files) == 0) {
  stop("No Stage5.csv files found under annotation_overlap/ — check the path is correct ",
       "and that you're running this from the same directory used for the annotation run.")
}

# xmsannotator's multilevelannotation() setwd()s into each output dir without
# restoring it, so later platform/db combos nest inside earlier ones' dirs
# (e.g. ".../c18/HMDB/annotation_overlap/hilic/KEGG/Stage5.csv"). The true
# platform/db combo is always the last two path components before Stage5.csv,
# not a prefix-stripped version of the full (nested) path.
# names like ".../c18/HMDB/annotation_overlap/hilic/KEGG/Stage5.csv" -> "hilic_KEGG"
stage5_names <- stage5_files %>%
  str_split("/") %>%
  map_chr(~ paste(.x[length(.x) - 2], .x[length(.x) - 1], sep = "_"))

stage5_tables <- stage5_files %>%
  map(read_csv, show_col_types = FALSE) %>%
  map(~ rename_all(.x, str_to_lower)) %>%
  set_names(stage5_names)

cat("\nStage5 table row counts:\n")
print(map_int(stage5_tables, nrow))

matched_all <- map2_dfr(stage5_tables, names(stage5_tables), function(tbl, nm) {
  mode <- str_extract(nm, "^[^_]+")   # "c18" or "hilic"
  lookup <- overlap_features %>% filter(platform == toupper(mode))
  map_dfr(seq_len(nrow(lookup)), function(i) {
    match_annotator(lookup$metabolite_feature[i], lookup$mz[i],
                     lookup$rt[i], tbl) %>%
      mutate(source = nm)
  }) %>% filter(match_chemical != "")
})

write_csv(matched_all, "overlap_metabolite_annotations_raw.csv")

best_match <- matched_all %>%
  group_by(id) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(metabolite_feature = id, annotated_compound = match_chemical)

write_csv(best_match, "overlap_metabolite_best_annotation.csv")

cat("\nFeatures annotated:", n_distinct(matched_all$id), "/", nrow(overlap_features), "\n")
cat("Written: overlap_metabolite_annotations_raw.csv and overlap_metabolite_best_annotation.csv\n")
