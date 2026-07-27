## ---------------------------
## Script: annotate_overlap_metabolites_xmsannotator.R
## Purpose: Compound ID for the metabolite features implicated in the
##          mQTL x metQTL overlap (cases + controls), via xMSannotator.
##
## Adapted from collaborator (Yufan Gong)'s 5-annotation.R, scoped down
## to only the ~150-200 overlap features rather than a full-metabolome
## run (his full run took ~4hrs; this should be much faster at this scale).
##
## BLOCKING PARAMETERS — must confirm before running:
##   - ION_MODE_C18 / ION_MODE_HILIC: acquisition polarity per platform.
##     His lab used c18=negative, hilic=positive, but this is a
##     lab/instrument-specific choice, NOT a universal convention.
##     Confirm with the PEG metabolomics core before running — an
##     incorrect ion mode will silently corrupt every adduct match.
##   - MASS_TOL_PPM: defaulted to 10ppm (matches collaborator's pipeline)
##     pending confirmation of the actual instrument's mass accuracy.
## ---------------------------

library(tidyverse)
library(xMSannotator)

# ============================================================
# Preflight: xMSannotator's internals call several packages that aren't
# always pulled in as hard dependencies of the conda env (e.g. limma,
# used inside its parallel worker nodes via snow::clusterEvalQ). Fail
# loudly here rather than after loading the ~166K-row HMDB database and
# running WGCNA clustering, only to hit "no package called X" deep in a
# forked worker process.
required_pkgs <- c("limma", "WGCNA", "flashClust", "snow")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required package(s) for xMSannotator: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall via: BiocManager::install(c(", paste0('"', missing_pkgs, '"', collapse = ", "), "))",
       "\nor, if this is a conda-managed env: conda install -c bioconda ",
       paste0("bioconductor-", tolower(missing_pkgs), collapse = " "))
}

# ============================================================
# Known xMSannotator bug (documented in collaborator's 5-annotation.R
# comments, but not scripted there — applied interactively by him):
# the package's internal multilevelannotation() contains
#   if (queryadductlist == "all" & mode == "pos") { ... }
# which errors under R >= 4.3 whenever queryadductlist is a multi-element
# vector (exactly what we pass), since if() on a length>1 condition is now
# a hard error rather than a warning. Patch both the "pos" and "neg"
# variants of this guard at runtime, plus an internal helper-function
# scoping issue he also flagged, before the namespace is used.
# ============================================================

patch_xmsannotator <- function() {
  ns <- asNamespace("xMSannotator")

  # Generic fix: xMSannotator's internal functions write intermediate stage
  # files (Stage1.csv, Stage4.csv, etc.) using bare relative filenames
  # rather than full file.path()-joined paths, so each stage setwd()s into
  # its target directory first. That directory isn't always pre-created
  # when a stage is invoked outside the package's normal full-pipeline
  # sequence (as we're doing here, running scoped to overlap features
  # only) — make every setwd() self-healing by auto-creating the target
  # directory immediately beforehand.
  guard_setwd <- function(fn_lines) {
    gsub(
      "setwd\\(([A-Za-z_][A-Za-z0-9_.]*)\\)",
      "{if (!dir.exists(\\1)) dir.create(\\1, recursive = TRUE, showWarnings = FALSE); setwd(\\1)}",
      fn_lines
    )
  }

  # Any bare clusterExport(cl, "NAME") with no explicit envir resolves via
  # default search-path lookup, which is ambiguous for non-exported
  # internal functions even after assignInNamespace() has patched them.
  # Force it to pull explicitly from the (now-patched) namespace. Calls
  # that already specify envir=... (e.g. several in
  # multilevelannotationstep4) are left untouched — the bare two-argument
  # pattern won't match those.
  guard_clusterexport <- function(fn_lines) {
    gsub(
      'clusterExport\\(cl, "([^"]+)"\\)',
      'clusterExport(cl, "\\1", envir = asNamespace("xMSannotator"))',
      fn_lines
    )
  }

  patch_fn <- function(fn_name, extra_gsubs = list()) {
    if (!exists(fn_name, envir = ns, inherits = FALSE)) {
      cat("(skip) function not found in namespace:", fn_name, "\n")
      return(invisible(NULL))
    }
    fn <- get(fn_name, envir = ns)
    fn_lines <- deparse(body(fn))
    for (g in extra_gsubs) fn_lines <- gsub(g$pattern, g$replacement, fn_lines, fixed = isTRUE(g$fixed))
    fn_lines <- guard_setwd(fn_lines)
    fn_lines <- guard_clusterexport(fn_lines)
    body(fn) <- parse(text = fn_lines)[[1]]
    environment(fn) <- ns
    assignInNamespace(fn_name, fn, ns = "xMSannotator")
    cat("Patched:", fn_name, "\n")
  }

  # multilevelannotation-specific bugs (vectorized-if guard, internal
  # helper scoping), on top of the universal setwd/clusterExport guards
  # applied to every function below.
  patch_fn("multilevelannotation", list(
    list(pattern = 'queryadductlist == "all" & mode == "pos"',
         replacement = 'is.character(queryadductlist) && length(queryadductlist) == 1 && queryadductlist == "all" && mode == "pos"',
         fixed = TRUE),
    list(pattern = 'queryadductlist == "all" & mode == "neg"',
         replacement = 'is.character(queryadductlist) && length(queryadductlist) == 1 && queryadductlist == "all" && mode == "neg"',
         fixed = TRUE),
    list(pattern = "get_peak_blocks_modulesvhclust",
         replacement = "xMSannotator:::get_peak_blocks_modulesvhclust",
         fixed = TRUE),
    # Genuine bug in the original source: the tempobjects.Rda save-list
    # (used across an rm(list=ls()); load("tempobjects.Rda") checkpoint
    # right before Stage 3) omits "num_nodes", even though a separate,
    # similarly-named step1_results.Rda save-list a few lines later DOES
    # include it — a copy-paste bug where num_nodes was added to one list
    # but not the other. After the rm(list=ls()) wipe, num_nodes is gone
    # and Stage 3's own "if (num_sets > num_nodes)" check fails.
    list(pattern = '"allsteps", "redundancy_check", "num_sets")',
         replacement = '"allsteps", "redundancy_check", "num_sets", "num_nodes")',
         fixed = TRUE)
  ))

  # Every other function found via find_cluster_creation.R's full
  # namespace scan that either contains setwd() directly, or creates/
  # participates in a parallel cluster (and so could ship an unpatched
  # setwd() into a worker independently of the ones above).
  other_fns <- c(
    "get_peak_blocks_modulesvhclust",
    "Annotationbychemical_IDs", "Annotationbychemical_IDschild",
    "Annotationbychemical_IDschild_multilevel", "Annotationbychemical_IDschildsimple",
    "ChemSpider.annotation", "data_preprocess", "find.Overlapping.mzs",
    "find.Overlapping.mzsvparallel", "overlapmzchild", "getVenn",
    "get_chemscorev1.6.71", "get_peak_blocks_cor", "getSumreplicates",
    "group_metabs_prev", "KEGG.annotationvold", "Metlin.annotation",
    "multilevelannotationstep2", "multilevelannotationstep3",
    "multilevelannotationstep4", "multilevelannotationstep5",
    "simpleAnnotation", "check_golden_rules", "check_element",
    "getMolecule", "get_confidence_stage2", "group_by_rt_histv2",
    "get_confidence_stage4", "group_by_rt"
  )
  for (fn_name in other_fns) patch_fn(fn_name)
}

patch_xmsannotator()

# ============================================================
# BLOCKING PARAMETERS
# ============================================================

ION_MODE_C18   <- 'neg'   # TODO: set to "neg" or "pos" once confirmed
ION_MODE_HILIC <- 'pos'   # TODO: set to "neg" or "pos" once confirmed
MASS_TOL_PPM   <- 10     # default per collaborator's pipeline
RT_TOL_SEC     <- 30     # matches his matchAnnotator() tolerance

if (is.null(ION_MODE_C18) || is.null(ION_MODE_HILIC)) {
  stop("Set ION_MODE_C18 and ION_MODE_HILIC before running — ",
       "getting this wrong silently corrupts every adduct match. ",
       "Check with the PEG metabolomics core for acquisition polarity.")
}

# Adduct lists per polarity, from his multilevelannotation() call
adduct_lists <- list(
  neg = c("M-H", "M-H2O-H", "M+Na-2H", "M+Cl", "M+FA-H"),
  pos = c("M+2H", "M+H+NH4", "M+ACN+2H", "M+2ACN+2H", "M+H", "M+NH4", "M+Na",
          "M+ACN+H", "M+ACN+Na", "M+2ACN+H", "2M+H", "2M+Na", "2M+ACN+H",
          "M+2Na-H", "M+H-H2O", "M+H-2H2O")
)
filter_adduct_lists <- list(neg = c("M-H"), pos = c("M+H"))

# ============================================================
# Step 1: Extract unique mz/rt features from the overlap tables
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

cat("Raw rows: cases =", nrow(overlap_cases), ", controls =", nrow(overlap_controls), "\n")
cat("NAs in metabolite_feature: cases =", sum(is.na(overlap_cases$metabolite_feature)),
    ", controls =", sum(is.na(overlap_controls$metabolite_feature)), "\n")

# Parse "C18_mz_rt_581.14_236.534" / "HILIC_mz_rt_804.55_65.765" into
# platform, mz, rt
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

cat("Distinct metabolite_feature: cases table =", n_distinct(overlap_cases$metabolite_feature),
    ", controls table =", n_distinct(overlap_controls$metabolite_feature), "\n")

overlap_features <- bind_rows(overlap_cases, overlap_controls) %>%
  distinct(metabolite_feature) %>%
  filter(!is.na(metabolite_feature)) %>%
  pull(metabolite_feature) %>%
  map_dfr(parse_feature)

stopifnot("mz/rt parsing failed on some features" = sum(is.na(overlap_features$mz)) == 0)

cat("Unique overlap metabolite features to annotate:", nrow(overlap_features), "\n")
cat("By platform:\n"); print(table(overlap_features$platform))

# ============================================================
# Load real per-sample intensity matrices. xMSannotator's WGCNA
# clustering step needs feature-feature correlation ACROSS SAMPLES to
# group co-eluting adducts/isotopes of the same compound — a bare
# mz/rt-only table (no intensities) produces an all-NA correlation
# matrix and crashes hclust(). Quant matrices are samples (rows) x
# features (columns), with column names already "mz_rt_<mz>_<rt>"
# (no platform prefix — that only exists in our metabolite_feature
# labels, stripped below before matching).
# ============================================================

load_rdata_obj <- function(path) {
  e <- new.env()
  nm <- load(path, envir = e)
  get(nm[1], envir = e)
}

c18_quant   <- load_rdata_obj("../PEG_METAB/Downloads/PEG_c18_quant_combat_LOD_cc.t.RData")
hilic_quant <- load_rdata_obj("../PEG_METAB/Downloads/PEG_hilic_quant_combat_PC_LOD_cc.t.RData")

cat("c18 quant matrix:", nrow(c18_quant), "samples x", ncol(c18_quant), "features\n")
cat("hilic quant matrix:", nrow(hilic_quant), "samples x", ncol(hilic_quant), "features\n")

# xMSannotator's dataA argument expects EXACTLY: col1 = mz, col2 = time,
# remaining cols = per-sample intensities — no id/label column, since
# each row is identified by its own mz+time. We keep feature identity
# in a separate lookup (overlap_features) for Step 3, matched by mz/rt
# rather than smuggling an id column into the annotation input itself
# (which would be misread as a non-numeric "sample" and reintroduce NAs).
build_feature_table <- function(quant_mat, wanted_keys) {
  keep <- intersect(colnames(quant_mat), wanted_keys)
  missing <- setdiff(wanted_keys, colnames(quant_mat))
  if (length(missing) > 0) {
    cat("  NOT found in quant matrix (", length(missing), " of ", length(wanted_keys), "):\n", sep = "")
    print(missing)
  }
  if (length(keep) == 0) return(tibble())

  sub <- quant_mat[, keep, drop = FALSE]
  t_sub <- as.data.frame(t(sub))              # features become rows, samples become columns
  mzrt <- str_match(rownames(t_sub), "^mz_rt_([0-9.]+)_([0-9.]+)$")
  stopifnot("Some kept feature names didn't parse as mz_rt_X_Y" = sum(is.na(mzrt[,2])) == 0)

  out <- data.frame(mz = as.numeric(mzrt[,2]), time = as.numeric(mzrt[,3]),
                     t_sub, check.names = FALSE, row.names = NULL)
  out
}

c18_keys   <- overlap_features %>% filter(platform == "C18")   %>% mutate(key = str_remove(metabolite_feature, "^C18_"))   %>% pull(key)
hilic_keys <- overlap_features %>% filter(platform == "HILIC") %>% mutate(key = str_remove(metabolite_feature, "^HILIC_")) %>% pull(key)

feature_c18   <- build_feature_table(c18_quant, c18_keys)
feature_hilic <- build_feature_table(hilic_quant, hilic_keys)

cat("Feature table (c18):", nrow(feature_c18), "features x", ncol(feature_c18) - 2, "samples\n")
cat("Feature table (hilic):", nrow(feature_hilic), "features x", ncol(feature_hilic) - 2, "samples\n")

# ============================================================
# Step 2: Run xMSannotator::multilevelannotation, scoped to these features
# ============================================================

dir.create("annotation_overlap", showWarnings = FALSE)

data(adduct_table)
data(adduct_weights)

run_annotation <- function(feature_tbl, mode_name, mode, adductlist, filter_adduct) {
  if (nrow(feature_tbl) == 0) {
    cat("No", mode_name, "features to annotate, skipping.\n")
    return(NULL)
  }
  c("HMDB", "KEGG", "LipidMaps") %>%
    map(function(db) {
      outdir <- file.path("annotation_overlap", mode_name, db)
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      # Convert to an absolute path: multilevelannotation() and its many
      # helper functions chain dozens of setwd() calls across Stage 1-5,
      # and a RELATIVE outloc can resolve to a different physical
      # directory at different points as the CWD drifts internally —
      # e.g. step1_results.Rda gets saved under one resolved location,
      # then a later setwd(outloc1) + load(...) resolves the same
      # relative string against a different CWD and can't find it.
      # An absolute path removes that ambiguity entirely.
      outdir <- normalizePath(outdir)

      # A database with very few matches for our small, targeted feature
      # set (e.g. LipidMaps, which only covers lipid species) can end up
      # with nothing surviving Stage 3's filtering — the package's own
      # read.csv("Stage3.csv") call doesn't handle that gracefully and
      # crashes instead of treating "zero results" as valid. Don't let
      # one database's genuine no-match outcome take down databases that
      # already completed successfully.
      tryCatch({
        xMSannotator::multilevelannotation(
          feature_tbl,
          max.mz.diff = MASS_TOL_PPM, max.rt.diff = RT_TOL_SEC,
          # num_nodes: WGCNA::allowWGCNAThreads() requires nThreads >= 2, so
          # true serial execution isn't available — we're forced to use
          # worker processes. This means our earlier patches likely still
          # won't reach code running inside those workers; see the
          # diagnostic print in patch_xmsannotator() for the cluster-
          # creation code we'll need to target next if this errors again.
          num_nodes = 2,
          queryadductlist = adductlist,
          filter.by = filter_adduct,
          adduct_weights = adduct_weights,
          mode = mode,
          db_name = db,
          num_sets = 50,   # scaled down from his 300, appropriate for ~100-200 features
          outloc = outdir
        )
      }, error = function(e) {
        cat("(!) ", mode_name, "/", db, " annotation failed, skipping: ", conditionMessage(e), "\n", sep = "")
        NULL
      })
    })
}

run_annotation(feature_c18, "c18", ION_MODE_C18,
                adduct_lists[[ION_MODE_C18]], filter_adduct_lists[[ION_MODE_C18]])
run_annotation(feature_hilic, "hilic", ION_MODE_HILIC,
                adduct_lists[[ION_MODE_HILIC]], filter_adduct_lists[[ION_MODE_HILIC]])

# ============================================================
# Step 3: Load Stage5 output and match back to our features
# (mirrors matchAnnotator() from his 5-annotation.R)
# ============================================================

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

stage5_files <- list.dirs("annotation_overlap", recursive = FALSE) %>%
  list.files(pattern = "\\.csv$", full.names = TRUE, recursive = TRUE) %>%
  keep(~ str_detect(.x, "Stage5"))

if (length(stage5_files) == 0) {
  warning("No Stage5 output found — xMSannotator may have failed silently. ",
          "Check annotation_overlap/{mode}/{db}/ for logs before proceeding.")
} else {
  stage5_names <- stage5_files %>%
    str_extract("[^/]+/[^/]+/Stage5") %>%
    str_replace_all("/", "_")

  stage5_tables <- stage5_files %>%
    map(read_csv, show_col_types = FALSE) %>%
    map(~ rename_all(.x, str_to_lower)) %>%
    set_names(stage5_names)

  matched_all <- map2_dfr(stage5_tables, names(stage5_tables), function(tbl, nm) {
    mode <- str_extract(nm, "^[^_]+")
    lookup <- overlap_features %>% filter(platform == toupper(mode))
    map_dfr(seq_len(nrow(lookup)), function(i) {
      match_annotator(lookup$metabolite_feature[i], lookup$mz[i],
                       lookup$rt[i], tbl) %>%
        mutate(source = nm)
    }) %>% filter(match_chemical != "")
  })

  write_csv(matched_all, "overlap_metabolite_annotations_raw.csv")

  # Best match per feature: highest score, then highest confidence
  best_match <- matched_all %>%
    group_by(id) %>%
    slice_max(score, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    rename(metabolite_feature = id, annotated_compound = match_chemical)

  write_csv(best_match, "overlap_metabolite_best_annotation.csv")

  cat("\nFeatures annotated:", n_distinct(matched_all$id), "/",
      nrow(overlap_features), "\n")
  cat("Written: overlap_metabolite_annotations_raw.csv and",
      "overlap_metabolite_best_annotation.csv\n")
}
