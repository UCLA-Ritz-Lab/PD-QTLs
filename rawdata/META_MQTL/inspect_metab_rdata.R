## Diagnostic: inspect structure of the metabolomics intensity matrices
## and the assay-ID <-> PEG-ID linking files, before building the
## feature table xMSannotator actually needs (with real per-sample
## intensities, not just mz/rt).

inspect_rdata <- function(path) {
  cat("\n==========", path, "==========\n")
  e <- new.env()
  loaded_names <- load(path, envir = e)
  cat("Objects loaded:", paste(loaded_names, collapse = ", "), "\n")

  for (nm in loaded_names) {
    obj <- get(nm, envir = e)
    cat("\n---", nm, "---\n")
    cat("Class:", paste(class(obj), collapse = "/"), "\n")
    if (is.matrix(obj) || is.data.frame(obj)) {
      cat("Dimensions:", paste(dim(obj), collapse = " x "), "\n")
      cat("First 5 row names:", paste(head(rownames(obj), 5), collapse = ", "), "\n")
      cat("First 5 col names:", paste(head(colnames(obj), 5), collapse = ", "), "\n")
      cat("Preview:\n")
      print(obj[1:min(5, nrow(obj)), 1:min(5, ncol(obj))])
    } else if (is.list(obj)) {
      cat("List length:", length(obj), "\n")
      cat("Names:", paste(head(names(obj), 10), collapse = ", "), "\n")
      str(obj, max.level = 1)
    } else {
      cat("Preview:\n")
      print(head(obj))
    }
  }
}

inspect_rdata("../PEG_METAB/Downloads/PEG_c18_quant_combat_LOD_cc.t.RData")
inspect_rdata("../PEG_METAB/Downloads/PEG_hilic_quant_combat_PC_LOD_cc.t.RData")
inspect_rdata("../PEG_METAB/Downloads/c18_keyVar_link_cc.RData")
inspect_rdata("../PEG_METAB/Downloads/hilic_keyVar_link_cc.RData")
