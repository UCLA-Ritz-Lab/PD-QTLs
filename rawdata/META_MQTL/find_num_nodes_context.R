library(xMSannotator)
ns <- asNamespace("xMSannotator")
fn_lines <- deparse(body(get("multilevelannotation", envir = ns)))

# Find every line mentioning num_nodes, and print with surrounding context
idx <- grep("num_nodes", fn_lines)
cat("Lines referencing num_nodes:", paste(idx, collapse = ", "), "\n\n")

for (i in idx) {
  lo <- max(1, i - 8)
  hi <- min(length(fn_lines), i + 8)
  cat("--- context around line", i, "---\n")
  cat(paste(fn_lines[lo:hi], collapse = "\n"), "\n\n")
}

# Also print context around num_sets / arg1,outloc1 closure specifically,
# since that's the parLapply call flagged in the earlier full-namespace scan
idx2 <- grep("num_sets|arg1, outloc1", fn_lines)
cat("\nLines referencing num_sets or the arg1/outloc1 closure:", paste(idx2, collapse=", "), "\n\n")
for (i in idx2) {
  lo <- max(1, i - 5)
  hi <- min(length(fn_lines), i + 20)
  cat("--- context around line", i, "---\n")
  cat(paste(fn_lines[lo:hi], collapse = "\n"), "\n\n")
}
