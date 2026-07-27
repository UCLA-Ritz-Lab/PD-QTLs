## Diagnostic: search the ENTIRE xMSannotator namespace for wherever the
## parallel cluster actually gets created (makeCluster/parLapply/etc.),
## since it isn't in multilevelannotation or get_peak_blocks_modulesvhclust
## themselves — the setwd() failure happens inside a worker process, so
## the cluster-creation call must live in some other internal function
## that calls into one of those (or is called by them).

library(xMSannotator)
ns <- asNamespace("xMSannotator")
all_fns <- ls(ns, all.names = TRUE)

cat("Total objects in namespace:", length(all_fns), "\n\n")

hits <- list()
for (fn_name in all_fns) {
  obj <- tryCatch(get(fn_name, envir = ns), error = function(e) NULL)
  if (!is.function(obj)) next
  fn_lines <- tryCatch(deparse(body(obj)), error = function(e) character(0))
  if (length(fn_lines) == 0) next

  cluster_hits <- grep("makeCluster|parLapply|clusterApply|clusterEvalQ|clusterExport|clusterCall|stopCluster",
                        fn_lines, value = TRUE)
  setwd_hits <- grep("setwd\\(", fn_lines, value = TRUE)

  if (length(cluster_hits) > 0) {
    cat("=== ", fn_name, " (cluster-related) ===\n")
    cat(paste(cluster_hits, collapse = "\n"), "\n\n")
    hits[[fn_name]] <- list(cluster = cluster_hits, setwd = setwd_hits)
  }
  if (length(setwd_hits) > 0 && length(cluster_hits) == 0) {
    cat("=== ", fn_name, " (setwd, no cluster call) ===\n")
    cat(paste(setwd_hits, collapse = "\n"), "\n\n")
  }
}

cat("\nFunctions containing BOTH a cluster call and a setwd():\n")
for (nm in names(hits)) {
  if (length(hits[[nm]]$setwd) > 0) cat(" -", nm, "\n")
}
