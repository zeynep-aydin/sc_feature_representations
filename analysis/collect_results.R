suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = "Collect all metadata files into a flat CSV")
parser$add_argument("--experiments_dir", default = NULL, help = "Path to experiments/ (default: PROJ_ROOT/experiments)")
parser$add_argument("--output", default = NULL, help = "Output CSV path (default: PROJ_ROOT/results/all_results_dev.csv)")
args <- parser$parse_args()

project_root <- Sys.getenv("PROJ_ROOT", unset = ".")
experiments_dir <- if (!is.null(args$experiments_dir)) args$experiments_dir else file.path(project_root, "experiments")
output_path <- if (!is.null(args$output)) args$output else file.path(project_root, "results", "all_results_dev.csv")

files <- list.files(experiments_dir, pattern = "_metadata[.]RData$", recursive = TRUE, full.names = TRUE)
cat(sprintf("Found %d metadata files in %s\n", length(files), experiments_dir))

.get <- function(x, default = NA) if (is.null(x)) default else x

extract_row <- function(f) {
  e <- new.env()
  load(f, envir = e)
  m <- e$metadata
  data.frame(
    file                  = basename(f),
    dataset               = .get(m$args$dataset),
    label_level           = .get(m$args$label_level),
    feature_mode          = .get(m$args$feature_mode, "intersection"),
    task                  = .get(m$args$task),
    method                = .get(m$method),
    algorithm             = .get(m$args$algorithm),
    train_pct             = .get(m$args$train_pct),
    n_dim                 = if (!is.null(m$args$n_dim)) m$args$n_dim else if (identical(m$method, "scimilarity")) .get(m$n_features) else NA,
    n_hvg                 = .get(m$args$n_hvg),
    max_epochs            = if (identical(.get(m$method), "scvi")) .get(m$args$max_epochs) else NA,
    run_id                = .get(m$run_id),
    n_train               = .get(m$n_train),
    n_test                = .get(m$n_test),
    actual_train_fraction = .get(m$actual_train_fraction),
    n_features            = .get(m$n_features),
    n_train_classes       = .get(m$n_train_classes),
    n_test_classes        = .get(m$n_test_classes),
    accuracy              = .get(m$metrics$accuracy),
    precision             = .get(m$metrics$precision),
    recall                = .get(m$metrics$recall),
    f1_macro              = .get(m$metrics$f1_score),
    f1_weighted           = .get(m$metrics$weighted_f1),
    auroc                 = .get(m$metrics$auroc),
    sigma                 = .get(m$rff_metadata$sigma),
    sigma_N               = .get(m$rff_metadata$sigma_N),
    preprocess_time_s     = .get(m$timing$preprocess_time),
    hvg_time_s            = .get(m$timing$hvg_time),
    reduction_time_s      = .get(m$timing$reduction_time),
    model_time_s          = .get(m$timing$model_training_time),
    total_time_s          = .get(m$timing$classification_script_time),
    peak_ram_gb           = .get(m$modeling_R_peak_gb),
    node_id               = .get(m$node_id),
    job_id                = .get(m$job_id),
    stringsAsFactors = FALSE
  )
}

rows <- lapply(files, function(f) {
  tryCatch(extract_row(f), error = function(e) {
    cat(sprintf("  SKIP %s: %s\n", basename(f), e$message))
    NULL
  })
})
results <- do.call(rbind, Filter(Negate(is.null), rows))

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(results, output_path, row.names = FALSE)
cat(sprintf("Written %d rows to %s\n", nrow(results), output_path))

# Summary table
cols <- c("dataset", "label_level", "task", "method", "n_dim", "n_hvg", "max_epochs", "run_id", "accuracy", "f1_macro", "auroc")
print(results[order(results$dataset, results$label_level, results$method, results$n_dim), intersect(cols, colnames(results))],
      digits = 4, row.names = FALSE)
