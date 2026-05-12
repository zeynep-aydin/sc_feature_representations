suppressPackageStartupMessages({
  library(Matrix)
  library(qs2)
})

DOUBLET_PATTERN <- paste(
  "doublet", "contaminant", "contamination", "low.quality", "low_quality", "unassigned", "unknown", "debris", "multiplet",
  sep = "|"
)

drop_by_label <- function(mtx, ct, don, extra_pattern = NULL) {
  pattern <- if (!is.null(extra_pattern)) paste(DOUBLET_PATTERN, extra_pattern, sep = "|") else DOUBLET_PATTERN
  keep <- !is.na(ct) & !grepl(pattern, ct, ignore.case = TRUE)
  cat(sprintf("Dropping %d / %d cells by label pattern\n", sum(!keep), length(keep)))
  list(mtx = mtx[keep, , drop = FALSE], ct = ct[keep], don = don[keep])
}

drop_by_min_cells <- function(mtx, ct, don, min_cells = 50) {
  counts <- table(ct)
  drop_types <- names(counts[counts < min_cells])
  if (length(drop_types) > 0) {
    drop_cells <- sum(counts[drop_types])
    cat(sprintf("Dropping %d type(s) with <%d cells (%d cells total): %s\n",
                length(drop_types), min_cells, drop_cells, paste(drop_types, collapse = ", ")))
  } else {
    cat(sprintf("No types dropped at min_cells=%d\n", min_cells))
  }
  keep <- !ct %in% drop_types
  list(mtx = mtx[keep, , drop = FALSE], ct = ct[keep], don = don[keep])
}

save_dataset <- function(mtx, ct, don, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  qs_save(mtx, file.path(out_dir, "counts.qs"))
  qs_save(as.factor(ct), file.path(out_dir, "labels.qs"))
  qs_save(as.character(don), file.path(out_dir, "donors.qs"))
  cat(sprintf("Saved %d cells x %d genes | %d types | %d donors -> %s\n",
              nrow(mtx), ncol(mtx),
              length(unique(ct)), length(unique(don)), out_dir))
}

export_h5ad <- function(mtx, out_dir, script_dir, conda_env = "scimilarity") {
  tmp_mtx   <- tempfile(fileext = ".mtx")
  tmp_genes <- tempfile(fileext = ".txt")
  out_h5ad  <- file.path(out_dir, "counts.h5ad")

  writeMM(mtx, tmp_mtx)
  writeLines(colnames(mtx), tmp_genes)

  cmd <- sprintf("conda run -n %s python %s --mtx %s --genes %s --out %s",
                 conda_env,
                 file.path(script_dir, "save_h5ad.py"),
                 tmp_mtx, tmp_genes, out_h5ad)
  ret <- system(cmd)
  unlink(c(tmp_mtx, tmp_genes))
  if (ret != 0) stop(sprintf("h5ad export failed (conda env: %s)", conda_env))
  cat(sprintf("Exported h5ad -> %s\n", out_h5ad))
  invisible(out_h5ad)
}
