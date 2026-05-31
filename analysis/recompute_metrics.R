suppressPackageStartupMessages({
  library(caret)
  library(pROC)
  library(argparse)
})

parser <- ArgumentParser(description = "Compute metrics for metadata files saved with --skip_metrics")
parser$add_argument("--experiments_dir", default = NULL)
parser$add_argument("--dry_run", action = "store_true", default = FALSE,
                    help = "Report which files would be patched without writing")
args <- parser$parse_args()

project_root <- Sys.getenv("PROJ_ROOT", unset = ".")
experiments_dir <- if (!is.null(args$experiments_dir)) args$experiments_dir else file.path(project_root, "experiments")

src_dir <- file.path(project_root, "src")
source(file.path(src_dir, "metrics.R"))

files <- list.files(experiments_dir, pattern = "_metadata[.]RData$", recursive = TRUE, full.names = TRUE)
cat(sprintf("Found %d metadata files in %s\n", length(files), experiments_dir))

n_patched <- 0
n_skipped <- 0
for (f in files) {
  e <- new.env()
  load(f, envir = e)
  m <- e$metadata
  if (!is.null(m$metrics)) {
    n_skipped <- n_skipped + 1
    next
  }
  cat(sprintf("  Patching: %s\n", basename(f)))
  if (!args$dry_run) {
    metrics <- tryCatch(
      compute_metrics(m$predicted_classes, m$class_probabilities, m$true_labels),
      error = function(err) { cat(sprintf("    ERROR: %s\n", err$message)); NULL }
    )
    if (!is.null(metrics)) {
      metadata <- m
      metadata$metrics <- metrics
      save(metadata, file = f)
      n_patched <- n_patched + 1
    }
  } else {
    n_patched <- n_patched + 1
  }
}
cat(sprintf("Patched %d files, skipped %d (already have metrics)\n", n_patched, n_skipped))
