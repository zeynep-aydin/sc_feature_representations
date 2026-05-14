suppressPackageStartupMessages({
  library(Matrix)
  library(MatrixExtra)
  library(sparseMatrixStats)
  library(qs2)
  library(rhdf5)
  library(argparse)
  library(RSpectra)
  library(glmnet)
  library(e1071)
  library(caret)
  library(pROC)
})

parser <- ArgumentParser(description = "Classification pipeline")
parser$add_argument("-r", "--run_id", type = "integer", help = "Run/replicate ID (seeds + split selection)")
parser$add_argument("-d", "--dataset", choices = c("zilionis_lung", "pbmc"))
parser$add_argument("-s", "--train_pct", default = 80L, type = "integer", choices = c(50L, 70L, 80L))
parser$add_argument("-m", "--method", choices = c("rff_lapl", "rff_gauss", "pca", "scimilarity", "scvi"))
parser$add_argument("-t", "--task", choices = c("tissue", "celltype"))
parser$add_argument("-a", "--algorithm", default = "glmnet", choices = c("glmnet", "svm"))
parser$add_argument("-n", "--n_dim", type = "integer", help = "Final output dimensions (RFF uses n_dim/2 internal projections)")
parser$add_argument("--n_hvg", type = "integer", help = "Top HVGs to select (optional)")
parser$add_argument("--max_epochs", type = "integer", help = "scVI max training epochs (default: 400)")
parser$add_argument("--label_level", choices = c("l1", "l2", "l3"), help = "Label granularity for multi-level datasets (e.g. pbmc)")
args <- parser$parse_args()

method <- if (!is.null(args$method)) args$method else "baseline"
run_id <- args$run_id
train_frac <- args$train_pct / 100
project_root <- Sys.getenv("PROJ_ROOT", unset = ".")
data_dir <- file.path(project_root, "data", args$dataset)

script_dir <- tryCatch(
  normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))),
  error = function(e) "dev/src"
)

source(file.path(script_dir, "utils.R"))
source(file.path(script_dir, "bandwidth.R"))
source(file.path(script_dir, "io.R"))
source(file.path(script_dir, "reduction.R"))
source(file.path(script_dir, "models.R"))
source(file.path(script_dir, "metrics.R"))

dataset_suffix <- if (!is.null(args$label_level)) {
  paste0(args$dataset, "_", args$label_level)
} else {
  args$dataset
}
split_label <- paste0("train", args$train_pct, "_run", run_id)
output_dir <- file.path(project_root, "experiments", dataset_suffix, method)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

########## PART 1: Load, normalize, HVG, scale ##########
log_info("Loading dataset...")
preprocess_start <- Sys.time()

raw <- load_dataset(args$dataset, data_dir, args$task, args$label_level)
split_info <- make_split(raw$mtx, raw$labels, raw$batch_labels, train_frac, run_id)

X_train <- raw$mtx[split_info$train_indices, , drop = FALSE]
X_test <- raw$mtx[split_info$test_indices, , drop = FALSE]
all_labels <- as.factor(raw$labels)
y_train <- all_labels[split_info$train_indices]
y_test <- all_labels[split_info$test_indices]
rm(raw)

missing_train <- setdiff(levels(y_train), unique(y_train))
missing_test <- setdiff(levels(y_test), unique(y_test))
if (length(missing_train) > 0)
  stop(sprintf("Train set missing classes: %s", paste(missing_train, collapse = ", ")))
if (length(missing_test) > 0)
  stop(sprintf("Test set missing classes: %s", paste(missing_test, collapse = ", ")))

hvg_time <- NULL
hvg_indices <- NULL

if (method != "scimilarity" && method != "scvi") {
  log_info("Log-normalizing datasets...")
  X_train <- log_normalize(X_train)
  X_test <- log_normalize(X_test)

  if (!is.null(args$n_hvg)) {
    hvg_start <- Sys.time()
    hvg_result <- select_hvg(X_train, X_test, args$n_hvg)
    X_train <- hvg_result$X_train
    X_test <- hvg_result$X_test
    hvg_indices <- hvg_result$hvg_indices
    hvg_time <- as.numeric(difftime(Sys.time(), hvg_start, units = "secs"))
    log_info(sprintf("HVG selection: %d genes retained (%.1fs)", args$n_hvg, hvg_time))
  }

  log_info("Min-max scaling datasets...")
  train_scaled <- min_max_scale(X_train)
  X_train <- train_scaled$X
  test_scaled <- min_max_scale(X_test, train_scaled$col_mins, train_scaled$col_maxs)
  X_test <- test_scaled$X
  rm(train_scaled, test_scaled)
}

preprocess_time <- as.numeric(difftime(Sys.time(), preprocess_start, units = "secs"))

########## PART 2: Dimensionality reduction ##########
reduction_seed <- NULL
rff_metadata <- NULL
filename_base <- "baseline"
reduction_time <- 0

if (method != "baseline") {
  result <- if (method == "rff_lapl") {
    rff_reduce(X_train, X_test, args$n_dim, run_id, kernel = "laplacian")
  } else if (method == "rff_gauss") {
    rff_reduce(X_train, X_test, args$n_dim, run_id, kernel = "gaussian")
  } else if (method == "pca") {
    pca_reduce(X_train, X_test, args$n_dim, run_id)
  } else if (method == "scvi") {
    scvi_embed(split_info, data_dir, project_root, run_id, args$n_dim, args$max_epochs)
  } else {
    scimilarity_embed(split_info, data_dir, project_root, run_id)
  }

  if (is.null(result$X_train) || is.null(result$X_test) ||
      ncol(result$X_train) != result$expected_dims || ncol(result$X_test) != result$expected_dims) {
    log_info(sprintf("ERROR: reduction failed! Expected %d dims, got train=%d test=%d",
                     result$expected_dims, ncol(result$X_train), ncol(result$X_test)))
    quit(status = 1)
  }

  X_train <- result$X_train
  X_test <- result$X_test
  rff_metadata <- result$rff_metadata
  filename_base <- result$filename_base
  reduction_seed <- result$reduction_seed
  reduction_time <- result$reduction_time
  if (!is.null(result$preprocess_time)) preprocess_time <- result$preprocess_time

  log_info(sprintf("Reduction successful: %d x %d", nrow(X_train), ncol(X_train)))
}

if (!is.null(args$n_hvg)) filename_base <- paste0("hvg", args$n_hvg, "_", filename_base)

########## PART 3: Train and evaluate ##########
log_info(sprintf("Training %s model...", args$algorithm))
model_seed <- 65 + run_id
set.seed(model_seed)
gc(reset = TRUE)

fit <- if (args$algorithm == "glmnet") train_glmnet(X_train, y_train) else train_svm(X_train, y_train)

preds <- if (args$algorithm == "glmnet") {
  predict_glmnet(fit$model, fit$best_lambda, X_test, y_test)
} else {
  predict_svm(fit$model, X_test, y_test)
}
metrics <- compute_metrics(preds$pred, preds$prob, y_test)

classification_script_time <- as.numeric(difftime(Sys.time(), preprocess_start, units = "secs"))
modeling_R_peak_gb <- get_peak_ram_gb()

########## PART 4: Save ##########
metadata <- list(
  node_id = Sys.info()["nodename"],
  job_id = Sys.getenv("SLURM_JOB_ID"),
  args = args,
  method = method,
  split_info = split_info,
  scaling_method = "min_max",
  n_train = nrow(X_train),
  n_test = nrow(X_test),
  actual_train_fraction = split_info$actual_train_fraction,
  run_id = run_id,
  n_features = ncol(X_test),
  n_train_classes = length(levels(y_train)),
  n_test_classes = length(levels(y_test)),
  train_sparsity = as.list(calculate_sparsity(X_train)),
  test_sparsity = as.list(calculate_sparsity(X_test)),
  reduction_seed = reduction_seed,
  model_seed = model_seed,
  rff_metadata = rff_metadata,
  hvg_indices = hvg_indices,
  metrics = metrics,
  predicted_classes = preds$pred,
  class_probabilities = preds$prob,
  true_labels = y_test,
  timing = list(
    preprocess_time = preprocess_time,
    hvg_time = hvg_time,
    reduction_time = reduction_time,
    model_training_time = fit$model_time,
    predict_time = preds$predict_time,
    classification_script_time = classification_script_time
  ),
  modeling_R_peak_gb = round(modeling_R_peak_gb, 3),
  filename_base = filename_base
)

model_object <- list(
  trained_model = fit$model,
  model_params = fit$model_params,
  training_context = list(
    r_version = R.version.string,
    package_versions = list(
      glmnet = if (args$algorithm == "glmnet") as.character(packageVersion("glmnet")) else NULL,
      e1071 = if (args$algorithm == "svm") as.character(packageVersion("e1071")) else NULL
    )
  )
)

metadata_filename <- paste0(split_label, "_", filename_base, "_", args$algorithm, "_", args$task, "_metadata.RData")
model_filename <- paste0(split_label, "_", filename_base, "_", args$algorithm, "_", args$task, "_model.RData")

save(metadata, file = file.path(output_dir, metadata_filename))
save(model_object, file = file.path(output_dir, model_filename))

log_info("Script completed successfully!")
