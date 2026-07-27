library(dplyr)
library(readr)
library(httr)
library(jsonlite)
library(purrr)

# ============================================================
# Annotate the mQTL x metQTL overlap SNPs with gene context via
# the Ensembl REST API. Metabolite features are left as their
# original mz_rt labels — collaborator ID mapping table pending.
# ============================================================

overlap_cases    <- read_csv("mqtl_metqtl_overlap_cases.csv",
                              col_types = cols(pair_id = col_character(), snp_id = col_character(),
                                                cpg_id = col_character(), status = col_character(),
                                                warnings = col_character(), category = col_character(),
                                                metabolite_feature = col_character(),
                                                same_sign_as_meta = col_logical(), replicated = col_logical(),
                                                dynamic = col_logical(), .default = col_double()),
                              guess_max = Inf)
overlap_controls <- read_csv("mqtl_metqtl_overlap_controls.csv",
                              col_types = cols(pair_id = col_character(), snp_id = col_character(),
                                                cpg_id = col_character(), status = col_character(),
                                                warnings = col_character(), category = col_character(),
                                                metabolite_feature = col_character(),
                                                same_sign_as_meta = col_logical(), replicated = col_logical(),
                                                dynamic = col_logical(), .default = col_double()),
                              guess_max = Inf)

stopifnot("snp_id parsing failed - NAs present in overlap_cases" = sum(is.na(overlap_cases$snp_id)) == 0)
stopifnot("snp_id parsing failed - NAs present in overlap_controls" = sum(is.na(overlap_controls$snp_id)) == 0)

unique_snps <- union(overlap_cases$snp_id, overlap_controls$snp_id)
cat("Unique SNPs to annotate:", length(unique_snps), "\n")

# ---- Ensembl REST lookup: overlapping gene, else nearest gene within window ----

WINDOW <- 250000   # bp, symmetric search window if no direct gene overlap
BASE_URL <- "https://rest.ensembl.org"

annotate_snp <- function(snp_id) {
  parts <- strsplit(snp_id, ":")[[1]]
  chr <- parts[1]
  pos <- as.integer(parts[2])

  query_region <- function(start, end) {
    url <- sprintf("%s/overlap/region/human/%s:%d-%d?feature=gene;content-type=application/json",
                    BASE_URL, chr, start, end)
    resp <- tryCatch(GET(url), error = function(e) NULL)
    if (is.null(resp) || status_code(resp) != 200) return(NULL)
    fromJSON(content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE)
  }

  # 1. Exact position first — is the SNP inside a gene body?
  hit <- query_region(pos, pos)
  Sys.sleep(0.1)   # stay well under Ensembl's rate limit

  if (!is.null(hit) && length(hit) > 0 && is.data.frame(hit) && nrow(hit) > 0) {
    gene_labels <- ifelse(is.na(hit$external_name) | hit$external_name == "",
                           hit$id, hit$external_name)
    genes <- unique(gene_labels)
    return(tibble(snp_id = snp_id, gene = paste(genes, collapse = ";"),
                   relationship = "overlapping", distance_bp = 0))
  }

  # 2. Widen the window and pick nearest by distance to gene start/end
  hit_wide <- query_region(max(1, pos - WINDOW), pos + WINDOW)
  Sys.sleep(0.1)

  if (is.null(hit_wide) || length(hit_wide) == 0 || !is.data.frame(hit_wide) || nrow(hit_wide) == 0) {
    return(tibble(snp_id = snp_id, gene = NA_character_,
                   relationship = paste0("none_within_", WINDOW, "bp"), distance_bp = NA_integer_))
  }

  hit_wide <- hit_wide %>%
    mutate(dist = pmax(0, pmax(start - pos, pos - end)),
           gene_label = ifelse(is.na(external_name) | external_name == "", id, external_name))

  nearest <- hit_wide %>% slice_min(dist, n = 1, with_ties = FALSE)

  tibble(snp_id = snp_id, gene = nearest$gene_label,
         relationship = "nearest", distance_bp = nearest$dist)
}

cat("Querying Ensembl REST API for", length(unique_snps), "SNPs (~", 
    round(length(unique_snps) * 0.2 / 60, 1), "min at current rate)...\n")

annotations <- map_dfr(unique_snps, function(s) {
  res <- tryCatch(annotate_snp(s), error = function(e) {
    tibble(snp_id = s, gene = NA_character_, relationship = "query_failed", distance_bp = NA_integer_)
  })
  res
})

cat("\nAnnotation outcome breakdown:\n")
print(table(annotations$relationship))

write_csv(annotations, "overlap_snp_gene_annotations.csv")
cat("Written: overlap_snp_gene_annotations.csv\n")

# ---- Join back onto the overlap result tables ----

overlap_cases_annot <- overlap_cases %>%
  left_join(annotations, by = "snp_id") %>%
  rename(nearest_gene = gene, gene_relationship = relationship, gene_distance_bp = distance_bp)

overlap_controls_annot <- overlap_controls %>%
  left_join(annotations, by = "snp_id") %>%
  rename(nearest_gene = gene, gene_relationship = relationship, gene_distance_bp = distance_bp)

write_csv(overlap_cases_annot, "mqtl_metqtl_overlap_cases_annotated.csv")
write_csv(overlap_controls_annot, "mqtl_metqtl_overlap_controls_annotated.csv")
cat("Written: mqtl_metqtl_overlap_cases_annotated.csv and mqtl_metqtl_overlap_controls_annotated.csv\n")

# ---- Quick summary: unique gene x SNP combos found ----
cat("\nUnique genes implicated (cases):\n")
print(overlap_cases_annot %>% distinct(snp_id, nearest_gene) %>% filter(!is.na(nearest_gene)) %>%
        count(nearest_gene, sort = TRUE))

cat("\nUnique genes implicated (controls):\n")
print(overlap_controls_annot %>% distinct(snp_id, nearest_gene) %>% filter(!is.na(nearest_gene)) %>%
        count(nearest_gene, sort = TRUE))

cat("\nNOTE: metabolite features remain as original mz_rt labels",
    "(e.g. C18_mz_rt_581.14_236.534) — pending collaborator compound ID mapping table.\n")
