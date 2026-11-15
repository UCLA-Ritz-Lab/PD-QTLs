#!/usr/bin/env Rscript
# ──────────────────────────────────────────────────────────────────────────────
# run_meta_mqtl.R — PEG + PPMI cis-mQTL meta-analysis (compute step)
#
# Split out of meta_mqtl_analysis.Rmd on 2026-11-14. After the PEG alignment and
# allele-coding fixes, the cohort files grew to ~40M rows at p < 1e-3 and the
# cross-cohort intersection to ~1.8M pairs, which no longer fits in a notebook
# on a workstation. This script does all the heavy work (meant for a
# high-memory cluster node) and writes compact outputs; meta_mqtl_analysis.Rmd
# now only reads them and draws tables/plots.
#
# Run from pipelines/META_MQTL:
#   Rscript run_meta_mqtl.R --cores 16
#
# Outputs (under --out_dir, default ../../results/META_MQTL):
#   meta/{group}_meta_results.csv.gz   full meta results per group, ranked by |z|, with
#                                      per-cohort beta_/se_ columns (+ in_peg2 for cases)
#   meta/{group}_sig_fdr05.csv         FDR_BH < 0.05 (read by PPMI_LONG; columns 1-3 stay
#                                      pair_id, snp_id, cpg_id — prep_full.sh reads 2-3 by position)
#   meta/{group}_lead_fdr05.csv        lead pair per CpG among FDR<0.05 (input to PPMI_LONG)
#   meta/{group}_top20.csv            top 20 pairs with per-study beta / p
#   qc/dataset_summary.csv             per-cohort counts and beta summary
#   qc/qq_points.csv.gz                thinned per-cohort QQ points
#   qc/meta_qq_points.csv.gz           thinned meta QQ points
#   qc/beta_sample.csv.gz              100k sampled betas per cohort
#   qc/run_info.txt                    settings, rows read, timings, sessionInfo
# ──────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(meta)        # metagen() for inverse-variance meta-analysis
  library(optparse)
  library(parallel)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--out_dir",      default = "../../results/META_MQTL"),
  make_option("--max_cis_p",    type = "double",  default = 1e-3,
              help = "common p threshold applied to every cohort on read [%default]"),
  make_option("--max_cis_rows", type = "double",  default = 15e6,
              help = "row-cap backstop per cohort file [%default]"),
  make_option("--cores",        type = "integer", default = 1L,
              help = "cores for the per-pair meta-analysis [%default]"),
  make_option("--engine",       default = "metagen",
              help = paste("'metagen' (one metagen() call per pair, ~0.06 s each) or",
                           "'closed_form' (exact vectorised common-effect IVW, Q and I2;",
                           "cross-checked against metagen() on 1000 pairs) [%default]")),
  make_option("--qq_grid",      type = "double",  default = 0.01,
              help = "QQ point thinning grid on the -log10 scale [%default]")
)))
stopifnot(opt$engine %in% c("metagen", "closed_form"))

t_start <- Sys.time()
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%F %T")), sprintf(...), "\n", sep = "")
META_DIR <- file.path(opt$out_dir, "meta"); QC_DIR <- file.path(opt$out_dir, "qc")
dir.create(META_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR,   recursive = TRUE, showWarnings = FALSE)

# ── Inputs and groupings ─────────────────────────────────────────────────────
RESULTS <- c(
  peg1_cases    = "../../results/PEG_MQTL/mqtl/peg1_cases/cis_mqtls.txt",
  peg1_controls = "../../results/PEG_MQTL/mqtl/peg1_controls/cis_mqtls.txt",
  peg2_cases    = "../../results/PEG_MQTL/mqtl/peg2_cases/cis_mqtls.txt",
  ppmi_cases    = "../../results/PPMI_MQTL/mqtl/ppmi_cases/cis_mqtls.txt",
  ppmi_controls = "../../results/PPMI_MQTL/mqtl/ppmi_controls/cis_mqtls.txt"
)
resolve_path <- function(p) if (file.exists(paste0(p, ".gz"))) paste0(p, ".gz") else p
RESULTS <- vapply(RESULTS, resolve_path, character(1))

# Primary: PEG1 x PPMI two-way, separately for cases and controls.
# Sensitivity: 3-way cases intersection including PEG2 (cases-only wave).
META_GROUPS        <- list(cases    = c("peg1_cases",    "ppmi_cases"),
                           controls = c("peg1_controls", "ppmi_controls"))
SENSITIVITY_GROUPS <- list(cases_peg2_sensitivity = c("peg1_cases", "peg2_cases", "ppmi_cases"))

# True number of cis tests per cohort, for QQ expected quantiles on truncated
# p-value sets (from MatrixEQTL SNP x CpG counts; July 2026, within ~0.02% of
# the re-run SNP counts).
n_cis_tests_df <- fread(file.path(QC_DIR, "n_cis_tests.csv"))
n_cis_tests    <- setNames(n_cis_tests_df$n_total, n_cis_tests_df$dataset)

# ── Load ─────────────────────────────────────────────────────────────────────
# PPMI_MQTL writes cis results at p < 1e-3, PEG_MQTL at p < 1e-2. Reading every
# cohort at the same threshold (default 1e-3) keeps them comparable and bounds
# memory; PEG pairs with 1e-3 <= p < 1e-2 do not enter the meta-analysis.
# MatrixEQTL output is p-sorted, so the row cap (a backstop) keeps the most
# significant rows. awk streams the gzip so no more than needed is expanded.
rows_read <- list()
load_mqtl <- function(path, nm) {
  log_msg("loading %s", nm)
  if (!file.exists(path)) stop("File not found: ", path)
  cmd <- sprintf("%s %s | awk -F'\\t' 'NR==1 || $5 < %g' | head -n %.0f",
                 if (grepl("\\.gz$", path)) "zcat" else "cat",
                 shQuote(path), opt$max_cis_p, opt$max_cis_rows + 1)
  d <- fread(cmd = cmd, colClasses = list(character = 1:2))
  if (nrow(d) == 0) stop(nm, ": read 0 rows — decompression or path problem?")
  setnames(d, c("snp_id", "cpg_id", "beta", "tstat", "pvalue", "fdr_dataset"))
  hit_cap <- nrow(d) == opt$max_cis_rows
  if (hit_cap) log_msg("*** %s HIT THE %.0f-ROW CAP (worst p %.3g) ***", nm, opt$max_cis_rows, max(d$pvalue))
  # PEG SNP IDs are chr:pos; older PPMI outputs were chr:pos:REF:ALT. Strip any
  # allele suffix so pair_id matches across cohorts. (Both pipelines now count
  # the ALT allele; see PEG_MQTL/PPMI_MQTL positions/export_allele.txt.)
  d[, snp_id := sub("^([0-9XYM]+:[0-9]+):[ACGT]+:[ACGT]+$", "\\1", snp_id)]
  d[, se := abs(beta / tstat)]
  d <- d[!is.na(se) & se > 0 & !is.na(beta) & !is.na(pvalue)]
  d[, pair_id := paste(snp_id, cpg_id, sep = "|")]
  d[, tstat := NULL]
  rows_read[[nm]] <<- data.table(dataset = nm, rows = nrow(d), hit_cap = hit_cap)
  log_msg("  %s: %s rows at p < %g", nm, format(nrow(d), big.mark = ","), opt$max_cis_p)
  d
}
mqtl_list <- lapply(names(RESULTS), function(nm) load_mqtl(RESULTS[[nm]], nm))
names(mqtl_list) <- names(RESULTS)

# ── Per-cohort QC outputs ────────────────────────────────────────────────────
log_msg("writing per-cohort QC")
fwrite(rbindlist(lapply(names(mqtl_list), function(nm) {
  d <- mqtl_list[[nm]]
  data.table(dataset = nm, n_assoc = nrow(d),
             n_fdr05 = sum(d$fdr_dataset < 0.05, na.rm = TRUE),
             n_p1e5 = sum(d$pvalue < 1e-5), n_p1e8 = sum(d$pvalue < 1e-8),
             median_beta = median(abs(d$beta)), max_beta = max(abs(d$beta)),
             n_extreme = sum(abs(d$beta) > 1))
})), file.path(QC_DIR, "dataset_summary.csv"))

# QQ: expected quantiles from ALL retained p-values (the smallest n of n_total
# tests, so ranks are exact); only the stored points are thinned, one per
# qq_grid cell on the -log10 scale.
# obs is -log10(p); pass neglog10p directly when p underflows (meta p from z).
qq_points <- function(p = NULL, n_total = NULL, grid = opt$qq_grid, neglog10p = NULL) {
  if (is.null(neglog10p)) {
    # MatrixEQTL floors p at ~2.2e-308; floor exact zeros too so ranks are kept
    p <- p[!is.na(p) & p <= 1]
    neglog10p <- -log10(pmax(p, .Machine$double.xmin))
  }
  obs <- sort(neglog10p[is.finite(neglog10p)], decreasing = TRUE)
  n <- length(obs); if (is.null(n_total)) n_total <- n
  r <- seq_len(n); a <- 0.5
  dt <- data.table(exp   = -log10((r - a) / (n_total + 1 - 2 * a)),
                   obs   = obs,
                   ci_hi = -log10(qbeta(0.025, r, n_total - r + 1)),
                   ci_lo = -log10(qbeta(0.975, r, n_total - r + 1)))
  dt <- dt[!duplicated(data.table(round(exp / grid), round(obs / grid)))]
  dt[, `:=`(n = n, n_total = n_total)]
  dt
}
fwrite(rbindlist(lapply(names(mqtl_list), function(nm)
  qq_points(mqtl_list[[nm]]$pvalue, n_cis_tests[[nm]])[, dataset := nm])),
  file.path(QC_DIR, "qq_points.csv.gz"))

set.seed(1)
fwrite(rbindlist(lapply(names(mqtl_list), function(nm) {
  d <- mqtl_list[[nm]]
  data.table(dataset = nm, beta = d$beta[sample.int(nrow(d), min(1e5, nrow(d)))])
})), file.path(QC_DIR, "beta_sample.csv.gz"))

# Exact -log10 two-sided p from z; meta p underflows to 0 below ~1e-308
neglog10p_z <- function(z) -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)

# ── Meta-analysis ────────────────────────────────────────────────────────────
# Common-effect (fixed) inverse-variance meta-analysis: cohorts are drawn from
# the same population and we want one common effect. SE = |beta / t|.
# Only pairs present in ALL constituent cohorts are meta-analysed.
meta_metagen <- function(B, S, studlab) {
  # one metagen() call per pair, chunked across cores; returns numeric columns
  idx    <- seq_len(nrow(B))
  chunks <- split(idx, ceiling(idx / 5000))
  one <- function(ii) t(vapply(ii, function(i) {
    mg <- tryCatch(metagen(TE = B[i, ], seTE = S[i, ], studlab = studlab,
                           method.tau = "DL", common = TRUE, random = FALSE),
                   error = function(e) NULL)
    if (is.null(mg)) return(rep(NA_real_, 7))
    c(mg$TE.common, mg$seTE.common, mg$pval.common, mg$zval.common,
      mg$Q, mg$pval.Q, mg$I2)
  }, numeric(7)))
  res <- if (opt$cores > 1) mclapply(chunks, one, mc.cores = opt$cores, mc.preschedule = FALSE)
         else lapply(chunks, one)
  bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(bad)) stop(sum(bad), " meta-analysis chunks failed (worker error or killed)")
  m <- do.call(rbind, res)
  colnames(m) <- c("meta_beta", "meta_se", "meta_pval", "meta_z", "Q", "Q_pval", "I2")
  as.data.table(m)
}
meta_closed_form <- function(B, S, studlab) {
  W <- 1 / S^2; k <- ncol(B)
  mb <- rowSums(W * B) / rowSums(W); ms <- sqrt(1 / rowSums(W)); z <- mb / ms
  Q  <- rowSums(W * (B - mb)^2)
  I2 <- pmax(0, (Q - (k - 1)) / Q); I2[!is.finite(I2)] <- 0
  out <- data.table(meta_beta = mb, meta_se = ms, meta_pval = 2 * pnorm(-abs(z)), meta_z = z,
                    Q = Q, Q_pval = pchisq(Q, k - 1, lower.tail = FALSE), I2 = I2)
  set.seed(1); chk <- sample(nrow(B), min(1000, nrow(B)))
  ref <- t(vapply(chk, function(i) {
    mg <- metagen(TE = B[i, ], seTE = S[i, ], studlab = studlab,
                  method.tau = "DL", common = TRUE, random = FALSE)
    c(mg$TE.common, mg$seTE.common, mg$pval.common, mg$Q, mg$pval.Q, mg$I2)
  }, numeric(6)))
  mine <- as.matrix(out[chk, .(meta_beta, meta_se, meta_pval, Q, Q_pval, I2)])
  stopifnot(isTRUE(all.equal(unname(mine), unname(ref), tolerance = 1e-6)))
  log_msg("  closed form matches metagen() on %d sampled pairs", length(chk))
  out
}

run_meta_group <- function(group, ds) {
  log_msg("=== %s: %s ===", group, paste(ds, collapse = " + "))
  common <- Reduce(intersect, lapply(mqtl_list[ds], `[[`, "pair_id"))
  log_msg("  pairs in all cohorts: %s", format(length(common), big.mark = ","))
  if (length(common) == 0) return(NULL)
  al <- lapply(ds, function(nm) {
    d <- mqtl_list[[nm]][pair_id %in% common, .(pair_id, snp_id, cpg_id, beta, se)]
    setnames(d, c("beta", "se"), paste0(c("beta_", "se_"), nm))
    if (nm != ds[1]) d[, c("snp_id", "cpg_id") := NULL]
    d
  })
  merged <- Reduce(function(a, b) merge(a, b, by = "pair_id"), al); rm(al)
  B <- as.matrix(merged[, paste0("beta_", ds), with = FALSE])
  S <- as.matrix(merged[, paste0("se_",   ds), with = FALSE])
  t0 <- Sys.time()
  res <- if (opt$engine == "metagen") meta_metagen(B, S, ds) else meta_closed_form(B, S, ds)
  log_msg("  meta (%s, %d cores) took %.1f min", opt$engine, opt$cores,
          as.numeric(Sys.time() - t0, units = "mins"))
  # keep per-cohort beta/se alongside the meta estimate
  out <- cbind(merged[, c("pair_id", "snp_id", "cpg_id", paste0(c("beta_", "se_"), rep(ds, each = 2))),
                      with = FALSE],
               res, n_datasets = length(ds))
  # rank by |z|, not p: meta p underflows to 0 for thousands of strong pairs,
  # which would leave their order arbitrary (same order as p otherwise)
  out <- out[!is.na(meta_pval)][order(-abs(meta_z))]
  out[, FDR_BH := p.adjust(meta_pval, method = "BH")]
  log_msg("  %s pairs; FDR<0.05: %s; p<1e-8: %s", format(nrow(out), big.mark = ","),
          format(sum(out$FDR_BH < 0.05), big.mark = ","), format(sum(out$meta_pval < 1e-8), big.mark = ","))
  out
}

meta_results <- lapply(names(META_GROUPS), function(g) run_meta_group(g, META_GROUPS[[g]]))
names(meta_results) <- names(META_GROUPS)
sens_results <- lapply(names(SENSITIVITY_GROUPS), function(g) run_meta_group(g, SENSITIVITY_GROUPS[[g]]))
names(sens_results) <- names(SENSITIVITY_GROUPS)

# Primary cases hits that are also present (p < max_cis_p) in PEG2
if (!is.null(meta_results$cases))
  meta_results$cases[, in_peg2 := pair_id %in% mqtl_list$peg2_cases$pair_id]

# ── Write meta outputs ───────────────────────────────────────────────────────
for (g in names(meta_results)) {
  d <- meta_results[[g]]; if (is.null(d)) next
  fwrite(d, file.path(META_DIR, paste0(g, "_meta_results.csv.gz")))
  # same format as the notebook wrote (write.csv, no quotes): read by PPMI_LONG
  write.csv(d[FDR_BH < 0.05], file.path(META_DIR, paste0(g, "_sig_fdr05.csv")),
            row.names = FALSE, quote = FALSE)
  # Lead pair per CpG (largest |z| among FDR<0.05): the input to PPMI_LONG.
  # Since the PEG fix nearly every intersected pair passes FDR (pairs enter only
  # if p < 1e-3 in every cohort), and SNPs in LD with the lead repeat its
  # signal, so one model per CpG instead of ~1.8M (decided 2026-11-14).
  # !duplicated keeps column order (prep_full.sh reads columns 2-3 by position).
  write.csv(d[FDR_BH < 0.05][order(-abs(meta_z))][!duplicated(cpg_id)],
            file.path(META_DIR, paste0(g, "_lead_fdr05.csv")),
            row.names = FALSE, quote = FALSE)
  # top 20 with per-study beta and p
  top <- d[1:min(20, .N)]
  for (nm in META_GROUPS[[g]]) {
    ps <- mqtl_list[[nm]][pair_id %in% top$pair_id, .(pair_id, beta, pvalue)]
    setnames(ps, c("beta", "pvalue"), paste0(nm, c("_beta", "_pvalue")))
    top <- merge(top, ps, by = "pair_id", all.x = TRUE, sort = FALSE)
  }
  fwrite(top[order(-abs(meta_z))], file.path(META_DIR, paste0(g, "_top20.csv")))
}
for (g in names(sens_results)) {
  d <- sens_results[[g]]; if (is.null(d)) next
  fwrite(d, file.path(META_DIR, paste0(g, "_meta_results.csv.gz")))
}
fwrite(rbindlist(lapply(names(meta_results), function(g)
  if (!is.null(meta_results[[g]]))
    qq_points(neglog10p = neglog10p_z(meta_results[[g]]$meta_z))[, group := g])),
  file.path(QC_DIR, "meta_qq_points.csv.gz"))

# ── Run info ─────────────────────────────────────────────────────────────────
sink(file.path(QC_DIR, "run_info.txt"))
cat("run_meta_mqtl.R finished", format(Sys.time()), "\n")
cat(sprintf("runtime: %.1f min\n", as.numeric(Sys.time() - t_start, units = "mins")))
cat("settings:\n"); str(opt[setdiff(names(opt), "help")])
cat("\nrows read per cohort:\n"); print(rbindlist(rows_read))
cat("\nmeta groups:\n")
for (g in c(names(meta_results), names(sens_results))) {
  d <- c(meta_results, sens_results)[[g]]
  cat(sprintf("  %-24s pairs=%s FDR<0.05=%s\n", g,
              if (is.null(d)) 0 else format(nrow(d), big.mark = ","),
              if (is.null(d)) 0 else format(sum(d$FDR_BH < 0.05), big.mark = ",")))
}
cat("\n"); print(sessionInfo())
sink()
log_msg("done in %.1f min", as.numeric(Sys.time() - t_start, units = "mins"))
