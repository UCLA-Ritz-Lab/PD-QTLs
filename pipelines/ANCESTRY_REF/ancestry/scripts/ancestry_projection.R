#!/usr/bin/env Rscript
# ──────────────────────────────────────────────────────────────────────────────
# Label projected ancestry PCs, export QTL covariates, and plot the
# 1000 Genomes / PPMI / PEG comparison.
#
# Every input here is a plink2 --score projection onto allele weights computed
# in the 1000 Genomes reference alone, so the reference and both cohorts share
# one coordinate system and can be plotted on common axes without rescaling.
#
# Self-reported race/ethnicity is used ONLY for the concordance report at the
# end. It never enters the PCA or the exported covariates.
#
# Usage (see ../Snakefile, rule ancestry_report):
#   Rscript ancestry_projection.R \
#     --kgp-sscore <f> --panel <f> --eigenval <f> --ref-eigenvec <f> \
#     --cohort NAME:sscore:selfreport_csv:id_column:id_prefix   (repeatable) \
#     --n-pcs-cov 10 --out-all <f> --out-covariate-dir <d> \
#     --plot-cohorts <png> --plot-facets <png> --plot-scree <png> \
#     --out-concordance <f> --out-projection-check <f>
# ──────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ── Palette ───────────────────────────────────────────────────────────────────
# Validated for scatter (all-pairs, light surface): worst pair CVD dE 24.7,
# normal-vision dE 33.6, both marks >= 3:1 against the surface.
#
# Only the two cohorts get a categorical hue. The 2504 reference samples are
# drawn in neutral grey and identified by direct centroid labels instead of a
# 5-hue legend -- five superpopulation hues plus two cohort hues would be seven
# simultaneously-visible categories in a scatter, which cannot clear the
# all-pairs colour-separation floors at any ordering.
PAL <- list(
  ppmi       = "#2a78d6",
  peg        = "#eb6834",
  reference  = "#4a3aa7",
  context    = "#c9c8c1",
  text_pri   = "#0b0b0b",
  text_sec   = "#52514e",
  grid       = "#e6e5e0",
  surface    = "#fcfcfb"
)
COHORT_COLOR <- c(PPMI = PAL$ppmi, PEG = PAL$peg)

theme_ancestry <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background   = element_rect(fill = PAL$surface, colour = NA),
      panel.background  = element_rect(fill = PAL$surface, colour = NA),
      panel.grid.major  = element_line(colour = PAL$grid, linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      axis.title        = element_text(colour = PAL$text_sec, size = base_size - 1),
      axis.text         = element_text(colour = PAL$text_sec, size = base_size - 2),
      plot.title        = element_text(colour = PAL$text_pri, face = "bold",
                                       size = base_size + 3),
      plot.subtitle     = element_text(colour = PAL$text_sec, size = base_size - 1),
      plot.caption      = element_text(colour = PAL$text_sec, size = base_size - 3,
                                       hjust = 0),
      legend.title      = element_text(colour = PAL$text_sec, size = base_size - 1),
      legend.text       = element_text(colour = PAL$text_pri, size = base_size - 1),
      legend.position   = "bottom",
      strip.text        = element_text(colour = PAL$text_pri, face = "bold",
                                       size = base_size - 1)
    )
}

# ── Argument parsing ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

get_opt <- function(flag, required = TRUE, default = NULL) {
  i <- which(args == flag)
  if (!length(i)) {
    if (required) stop("Missing required argument: ", flag, call. = FALSE)
    return(default)
  }
  args[i[1] + 1L]
}
get_all <- function(flag) {
  i <- which(args == flag)
  if (!length(i)) character(0) else args[i + 1L]
}

kgp_sscore_f <- get_opt("--kgp-sscore")
panel_f      <- get_opt("--panel")
eigenval_f   <- get_opt("--eigenval")
ref_eigvec_f <- get_opt("--ref-eigenvec")
n_pcs_cov    <- as.integer(get_opt("--n-pcs-cov", required = FALSE, default = "10"))
out_all_f    <- get_opt("--out-all")
out_cov_dir  <- get_opt("--out-covariate-dir")
plot_coh_f   <- get_opt("--plot-cohorts")
plot_fac_f   <- get_opt("--plot-facets")
plot_scr_f   <- get_opt("--plot-scree")
out_conc_f   <- get_opt("--out-concordance")
out_proj_f   <- get_opt("--out-projection-check")

cohort_specs <- lapply(get_all("--cohort"), function(s) {
  p <- strsplit(s, ":", fixed = TRUE)[[1]]
  if (length(p) != 5L)
    stop("--cohort must be NAME:sscore:selfreport:id_column:id_prefix, got: ", s,
         call. = FALSE)
  list(name = p[1], sscore = p[2], selfreport = p[3],
       id_col = p[4], id_prefix = if (p[5] == "NONE") "" else p[5])
})
if (!length(cohort_specs)) stop("No --cohort specifications given", call. = FALSE)

# Optional per-cohort ID crosswalk, as NAME:csv with columns vcf_id,cohort_id.
#
# A cohort needs one when its covariate files are keyed differently from its
# genotypes: PEG's .sscore rows are GWAS_IDs (CRG_*) while ancestry_peg.csv is
# keyed by Pegid, so the self-report join has nothing to match on and would
# quietly label the entire cohort "Unknown". Written by cohort_keep_list.
id_maps <- lapply(get_all("--id-map"), function(s) {
  p <- strsplit(s, ":", fixed = TRUE)[[1]]
  if (length(p) != 2L)
    stop("--id-map must be NAME:csv, got: ", s, call. = FALSE)
  p[2]
})
names(id_maps) <- vapply(get_all("--id-map"), function(s)
  strsplit(s, ":", fixed = TRUE)[[1]][1], "")

# ── Readers ───────────────────────────────────────────────────────────────────

#' Read a plink2 .sscore projection into IID + PC1..PCn.
#'
#' plink2 names score columns "<weightcol>_AVG"; because the weights file
#' carries PC1..PCn headers, the projections arrive as PC1_AVG..PCn_AVG. The
#' _AVG form is what we want: it is the mean over *non-missing* genotypes, so
#' samples with patchier coverage are not dragged toward the origin the way a
#' raw sum would drag them.
read_sscore <- function(path) {
  dt <- fread(path)
  setnames(dt, sub("^#", "", names(dt)))
  if (!"IID" %in% names(dt)) stop("No IID column in ", path, call. = FALSE)
  pc_cols <- grep("^PC[0-9]+_AVG$", names(dt), value = TRUE)
  if (!length(pc_cols)) pc_cols <- grep("^PC[0-9]+$", names(dt), value = TRUE)
  if (!length(pc_cols)) stop("No PC score columns found in ", path, call. = FALSE)
  # Order numerically (PC10 must not sort before PC2).
  pc_cols <- pc_cols[order(as.integer(sub("\\D*([0-9]+).*", "\\1", pc_cols)))]
  out <- dt[, c("IID", pc_cols), with = FALSE]
  setnames(out, c("IID", paste0("PC", seq_along(pc_cols))))
  out[, IID := as.character(IID)]
  out[]
}

read_eigenvec <- function(path) {
  dt <- fread(path)
  setnames(dt, sub("^#", "", names(dt)))
  pc_cols <- grep("^PC[0-9]+$", names(dt), value = TRUE)
  pc_cols <- pc_cols[order(as.integer(sub("\\D*([0-9]+).*", "\\1", pc_cols)))]
  out <- dt[, c("IID", pc_cols), with = FALSE]
  out[, IID := as.character(IID)]
  out[]
}

#' Collapse a cohort's self-reported race/ethnicity to one label.
#'
#' Detected from the file's own columns rather than hardcoded per cohort:
#'   * PEG  supplies a single categorical column (`Race_gwas`).
#'   * PPMI supplies PPMI-style multi-select binary indicators (RAWHITE,
#'     RABLACK, RAASIAN, RAINDALS, RAHAWOPI, RANOS, RAUNKNOWN) plus HISPLAT.
#'     Multi-select means a participant can tick several boxes, so anyone with
#'     more than one flag collapses to "More than one race" rather than being
#'     silently assigned to whichever column happens to be checked first.
read_self_report <- function(path, id_col, id_prefix) {
  dt <- fread(path)
  if (!id_col %in% names(dt))
    stop("id column '", id_col, "' not in ", path, call. = FALSE)
  dt[, .join_id := paste0(id_prefix, trimws(as.character(get(id_col))))]

  ra_map <- c(RAWHITE = "White", RABLACK = "Black", RAASIAN = "Asian",
              RAINDALS = "American Indian/Alaska Native",
              RAHAWOPI = "Native Hawaiian/Pacific Islander", RANOS = "Other")
  ra_present <- intersect(names(ra_map), names(dt))

  if ("Race_gwas" %in% names(dt)) {
    lab <- as.character(dt$Race_gwas)
    lab[is.na(lab) | lab == ""] <- "Unknown"
  } else if (length(ra_present) > 0L) {
    m <- as.matrix(dt[, lapply(.SD, function(x) as.integer(x == 1)),
                      .SDcols = ra_present])
    m[is.na(m)] <- 0L
    n_flag <- rowSums(m)
    lab <- ifelse(
      n_flag > 1L, "More than one race",
      ifelse(n_flag == 1L, ra_map[ra_present][max.col(m, ties.method = "first")],
             "Unknown"))
    if ("RAUNKNOWN" %in% names(dt))
      lab[!is.na(dt$RAUNKNOWN) & dt$RAUNKNOWN == 1 & n_flag == 0L] <- "Unknown"
    # Ethnicity is a separate axis from race in this coding. Only 1 is
    # documented as "Yes"; other non-zero codes are left out of the label
    # rather than guessed at.
    if ("HISPLAT" %in% names(dt)) {
      hisp <- !is.na(dt$HISPLAT) & dt$HISPLAT == 1
      lab[hisp] <- paste0(lab[hisp], " (Hispanic/Latino)")
    }
  } else {
    warning("No recognised race columns in ", path, "; labelling all Unknown")
    lab <- rep("Unknown", nrow(dt))
  }

  unique(data.table(IID = dt$.join_id, self_report = lab), by = "IID")
}

# ── Load ──────────────────────────────────────────────────────────────────────
message("Reading reference projection: ", kgp_sscore_f)
kgp <- read_sscore(kgp_sscore_f)

# header = TRUE is required, not cosmetic: the 1000G panel's header line ends in
# two trailing tabs, so it declares six columns while every data row has four.
# With fill = TRUE alone, fread cannot tell that first line is a header, keeps it
# as a data row, and names the columns V1..V6 — at which point the sample and
# super_pop lookups below fail. Forcing the header, then dropping the two
# phantom columns and the empty row, recovers the real 2504-sample table.
panel <- fread(panel_f, fill = TRUE, header = TRUE)
setnames(panel, tolower(names(panel)))
panel <- panel[, which(names(panel) != "" &
                       !grepl("^v[0-9]+$", names(panel))), with = FALSE]
stopifnot("panel needs sample/super_pop columns" =
            all(c("sample", "super_pop") %in% names(panel)))
panel[, sample := as.character(sample)]
panel <- panel[!is.na(sample) & sample != ""]

kgp <- merge(kgp, panel[, .(IID = sample, population = pop,
                            superpopulation = super_pop)],
             by = "IID", all.x = TRUE)
kgp[is.na(superpopulation), superpopulation := "Unknown"]
kgp[, `:=`(cohort = "1000G", self_report = NA_character_)]

cohort_dts <- lapply(cohort_specs, function(sp) {
  message("Reading cohort ", sp$name, ": ", sp$sscore)
  d <- read_sscore(sp$sscore)
  sr <- read_self_report(sp$selfreport, sp$id_col, sp$id_prefix)

  # Translate self-report IDs into the .sscore's ID space where the two differ.
  if (!is.null(id_maps[[sp$name]])) {
    im <- fread(id_maps[[sp$name]], colClasses = "character")
    for (needed in c("vcf_id", "cohort_id"))
      if (!needed %in% names(im))
        stop("id-map ", id_maps[[sp$name]], " lacks a ", needed, " column",
             call. = FALSE)
    sr <- merge(sr, im, by.x = "IID", by.y = "cohort_id")
    sr[, IID := vcf_id][, vcf_id := NULL]
    message(sprintf("  %s: id-map translated %d self-report IDs", sp$name, nrow(sr)))
  }

  d <- merge(d, sr, by = "IID", all.x = TRUE)
  d[is.na(self_report), self_report := "Unknown"]

  # An all-Unknown cohort means the join found nothing — almost always an ID
  # namespace mismatch, which is silent otherwise. Warn loudly.
  if (nrow(d) > 0L && all(d$self_report == "Unknown"))
    warning(sp$name, ": no self-report IDs matched the projection; ",
            "check the id_column / id-map for this cohort", call. = FALSE)
  d[, `:=`(cohort = sp$name, population = NA_character_,
           superpopulation = NA_character_)]
  d[]
})
names(cohort_dts) <- vapply(cohort_specs, `[[`, "", "name")

n_pc_common <- min(vapply(c(list(kgp), cohort_dts),
                          function(d) sum(grepl("^PC[0-9]+$", names(d))),
                          integer(1)))
pc_cols <- paste0("PC", seq_len(n_pc_common))
message("PCs available across all inputs: ", n_pc_common)

keep <- c("IID", "cohort", "population", "superpopulation", "self_report", pc_cols)
scores <- rbindlist(c(list(kgp[, ..keep]), lapply(cohort_dts, function(d) d[, ..keep])),
                    use.names = TRUE)

for (nm in names(cohort_dts))
  message(sprintf("  %-6s: %d samples", nm, nrow(cohort_dts[[nm]])))
message(sprintf("  %-6s: %d samples", "1000G", nrow(kgp)))

# ── Variance explained ────────────────────────────────────────────────────────
eigenval <- scan(eigenval_f, quiet = TRUE)
# plink2 emits only the requested top-N eigenvalues, so this is the share of
# variance among those N -- not of total genotypic variance. Labelled as such.
pct_var <- round(100 * eigenval / sum(eigenval), 1)
axlab <- function(i) sprintf("PC%d (%.1f%% of top-%d variance)",
                             i, pct_var[i], length(eigenval))

# ── Projection sanity check ───────────────────────────────────────────────────
# Projecting the reference through its own weights must reproduce the directly
# computed eigenvectors up to a per-PC scale (and arbitrary sign). Anything
# below |r| ~ 0.99 means the weights, the allele counts, or the variant set
# used for scoring disagree with the ones the PCA was built on.
ref_eig <- read_eigenvec(ref_eigvec_f)
chk <- merge(kgp[, c("IID", pc_cols), with = FALSE], ref_eig,
             by = "IID", suffixes = c(".proj", ".direct"))
proj_lines <- c(
  "Projection sanity check: reference projected through its own allele weights",
  "vs. the directly computed eigenvectors. |r| should be ~1.000 for every PC.",
  "(Sign is arbitrary; scale differs because --score reports a per-allele mean.)",
  "",
  sprintf("%-6s %10s %12s", "PC", "|r|", "status")
)
worst_r <- 1
for (i in seq_len(min(n_pc_common, ncol(ref_eig) - 1L))) {
  a <- chk[[paste0("PC", i, ".proj")]]; b <- chk[[paste0("PC", i, ".direct")]]
  if (is.null(a) || is.null(b)) next
  r <- abs(suppressWarnings(cor(a, b)))
  if (is.na(r)) r <- 0
  worst_r <- min(worst_r, r)
  proj_lines <- c(proj_lines,
                  sprintf("%-6s %10.4f %12s", paste0("PC", i), r,
                          if (r >= 0.99) "ok" else "CHECK"))
}
proj_lines <- c(proj_lines, "",
                sprintf("worst |r| across PCs: %.4f", worst_r),
                if (worst_r >= 0.99) "RESULT: PASS"
                else "RESULT: CHECK - projection does not reproduce the reference PCA")
writeLines(proj_lines, out_proj_f)
message(paste(tail(proj_lines, 2), collapse = " "))

# ── Nearest-superpopulation assignment ────────────────────────────────────────
# Assign each cohort sample to the nearest 1000G superpopulation centroid in
# the first few PCs. This is a coarse descriptive label for QC only -- it is a
# hard assignment over continuous axes, so admixed samples land at whichever
# centroid is nearest and genuinely intermediate ancestry is not represented.
# It exists to be cross-tabulated against self-report, not to be used as data.
n_assign_pc <- min(4L, n_pc_common)
apc <- paste0("PC", seq_len(n_assign_pc))
centroids <- kgp[superpopulation != "Unknown",
                 lapply(.SD, mean), by = superpopulation, .SDcols = apc]
cmat <- as.matrix(centroids[, ..apc])
rownames(cmat) <- centroids$superpopulation

nearest_superpop <- function(dt) {
  m <- as.matrix(dt[, ..apc])
  d <- sapply(seq_len(nrow(cmat)), function(k)
    rowSums((m - matrix(cmat[k, ], nrow(m), ncol(m), byrow = TRUE))^2))
  if (is.null(dim(d))) d <- matrix(d, nrow = 1)
  rownames(cmat)[max.col(-d, ties.method = "first")]
}
scores[, assigned_superpop := NA_character_]
for (nm in names(cohort_dts))
  scores[cohort == nm, assigned_superpop := nearest_superpop(.SD), .SDcols = apc]
scores[cohort == "1000G", assigned_superpop := superpopulation]

# ── Self-report vs PCA concordance ────────────────────────────────────────────
conc <- c(
  "Self-reported race/ethnicity vs. nearest 1000 Genomes superpopulation centroid",
  sprintf("(assignment uses PC1-PC%d; self-report never entered the PCA)", n_assign_pc),
  "",
  "Discordance here is expected and is not by itself an error: self-reported",
  "race is a social variable, the centroid assignment is a hard cut over",
  "continuous axes, and admixed individuals have no correct single answer.",
  "Read it for gross problems -- a whole cohort shifted, or a sample swap.",
  ""
)
for (nm in names(cohort_dts)) {
  d <- scores[cohort == nm]
  tb <- table(`self-reported` = d$self_report, `nearest 1000G` = d$assigned_superpop)
  conc <- c(conc, paste0("── ", nm, " (n = ", nrow(d), ") ",
                         strrep("─", max(1, 50 - nchar(nm)))), "",
            capture.output(print(tb)), "")
}
writeLines(conc, out_conc_f)

# ── Exports ───────────────────────────────────────────────────────────────────
fwrite(scores, out_all_f)
message("Wrote ", out_all_f, " (", nrow(scores), " samples)")

cov_cols <- paste0("PC", seq_len(min(n_pcs_cov, n_pc_common)))
for (nm in names(cohort_dts)) {
  f <- file.path(out_cov_dir, paste0("ancestry_pcs_", nm, ".csv"))
  fwrite(scores[cohort == nm, c("IID", cov_cols), with = FALSE], f)
  message("Wrote ", f, " (", nrow(scores[cohort == nm]), " samples, ",
          length(cov_cols), " PCs)")
}

# ── Plots ─────────────────────────────────────────────────────────────────────
scores[, plot_group := factor(cohort, levels = c("1000G", names(cohort_dts)))]
ref_pts <- scores[cohort == "1000G"]
coh_pts <- scores[cohort != "1000G"]

cent_lab <- ref_pts[superpopulation != "Unknown",
                    .(PC1 = median(PC1), PC2 = median(PC2),
                      PC3 = if ("PC3" %in% names(ref_pts)) median(PC3) else NA_real_,
                      PC4 = if ("PC4" %in% names(ref_pts)) median(PC4) else NA_real_),
                    by = superpopulation]

# ── Plot 1: cohorts against the reference backdrop ────────────────────────────
# The reference is context, so it is drawn in neutral grey and named by direct
# centroid labels; only the two cohorts carry a categorical hue.
#' @param label_centroids draw superpopulation names at their centroids.
#'   Only meaningful on PC1/PC2, where the superpopulations actually separate.
#'   On later PCs every centroid collapses toward the origin and the labels
#'   land on top of each other, so they are suppressed there rather than
#'   rendered as an unreadable pile.
panel_plot <- function(xv, yv, label_centroids = FALSE) {
  ggplot() +
    geom_point(data = ref_pts, aes(.data[[xv]], .data[[yv]]),
               colour = PAL$context, size = 0.7, alpha = 0.6) +
    geom_point(data = coh_pts, aes(.data[[xv]], .data[[yv]], colour = cohort),
               size = 1.1, alpha = 0.75) +
    {if (label_centroids)
      geom_label(data = cent_lab,
                 aes(.data[[xv]], .data[[yv]], label = superpopulation),
                 colour = PAL$text_sec, fill = PAL$surface, alpha = 0.72,
                 label.size = 0, label.padding = unit(0.12, "lines"),
                 size = 3.4, fontface = "bold")
     else NULL} +
    scale_colour_manual(values = COHORT_COLOR, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1))) +
    labs(x = axlab(as.integer(sub("PC", "", xv))),
         y = axlab(as.integer(sub("PC", "", yv)))) +
    theme_ancestry()
}

p1 <- panel_plot("PC1", "PC2", label_centroids = TRUE)
have34 <- all(c("PC3", "PC4") %in% names(scores))
p2 <- if (have34) panel_plot("PC3", "PC4") else NULL

png(plot_coh_f, width = 2200, height = 1150, res = 170)
if (!is.null(p2) && requireNamespace("patchwork", quietly = TRUE)) {
  print(
    patchwork::wrap_plots(p1, p2, ncol = 2, guides = "collect") +
      patchwork::plot_annotation(
        title = "PPMI and PEG projected onto 1000 Genomes ancestry PCs",
        subtitle = paste0("Grey points are the ", nrow(ref_pts),
                          " 1000 Genomes phase 3 reference samples; labels mark ",
                          "superpopulation centroids on the PC1/PC2 panel"),
        caption = paste("PCs computed in the reference alone and applied to both",
                        "cohorts as fixed allele weights, so all three sets share",
                        "one coordinate system.\nGRCh37/hg19."),
        theme = theme_ancestry()
      ) & theme(legend.position = "bottom")
  )
} else {
  print(p1 + labs(title = "PPMI and PEG projected onto 1000 Genomes ancestry PCs"))
}
invisible(dev.off())
message("Wrote ", plot_coh_f)

# ── Plot 2: small multiples, one group per facet ──────────────────────────────
# Small multiples rather than seven hues in one scatter: each facet highlights a
# single group against the full reference in grey, so identity comes from the
# strip label and no two categorical hues ever compete inside one panel.
facet_groups <- c(sort(unique(ref_pts[superpopulation != "Unknown"]$superpopulation)),
                  names(cohort_dts))
fac <- rbindlist(lapply(facet_groups, function(g) {
  d <- if (g %in% names(cohort_dts)) scores[cohort == g]
       else ref_pts[superpopulation == g]
  if (!nrow(d)) return(NULL)
  cbind(d[, c("PC1", "PC2"), with = FALSE], facet_group = g)
}))
fac[, facet_group := factor(facet_group, levels = facet_groups)]
fac[, hue := fifelse(facet_group %in% names(cohort_dts),
                     COHORT_COLOR[as.character(facet_group)], PAL$reference)]

backdrop <- scores[, .(PC1, PC2)]

p_fac <- ggplot() +
  geom_point(data = backdrop, aes(PC1, PC2), colour = PAL$context,
             size = 0.35, alpha = 0.5) +
  geom_point(data = fac, aes(PC1, PC2, colour = hue), size = 0.8, alpha = 0.85) +
  scale_colour_identity() +
  facet_wrap(~ facet_group, ncol = 4) +
  labs(
    title = "Where each group sits in the reference PC space",
    subtitle = "Each panel highlights one group; grey shows all samples for context",
    x = axlab(1), y = axlab(2),
    caption = "1000 Genomes phase 3 superpopulations (AFR, AMR, EAS, EUR, SAS) followed by the two study cohorts."
  ) +
  theme_ancestry()

png(plot_fac_f, width = 2200, height = 1250, res = 170)
print(p_fac)
invisible(dev.off())
message("Wrote ", plot_fac_f)

# ── Plot 3: scree ─────────────────────────────────────────────────────────────
scree <- data.table(PC = seq_along(pct_var), pct = pct_var)
scree[, lab := ifelse(PC <= 10, sprintf("%.1f", pct), NA_character_)]

p_scree <- ggplot(scree, aes(factor(PC), pct)) +
  geom_col(fill = PAL$ppmi, width = 0.68) +
  geom_text(aes(label = lab), vjust = -0.5, colour = PAL$text_sec,
            size = 2.9, na.rm = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Variance captured by each reference PC",
    subtitle = sprintf("Share of variance among the top %d PCs; first %d are exported as covariates",
                       length(pct_var), length(cov_cols)),
    x = "Principal component", y = "% of top-PC variance"
  ) +
  theme_ancestry() +
  theme(panel.grid.major.x = element_blank())

png(plot_scr_f, width = 1700, height = 950, res = 170)
print(p_scree)
invisible(dev.off())
message("Wrote ", plot_scr_f)

message("Done.")
