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

gc(reset = TRUE)

parser <- ArgumentParser(description = "Classification pipeline")

parser$add_argument("-r", "--replicate", default = 1, type = "integer", help = "Replicate number")
parser$add_argument("-d", "--dataset", choices = c("scea", "zhao_immune", "zilionis_lung", "he_organs"), help = "Dataset to use")
parser$add_argument("-f", "--feature_mode", default = "intersection", choices = c("intersection", "thresh08"), help = "(For SCEA dataset)")
parser$add_argument("-s", "--train_split", default = 0.8, type = "double", help = "Train split ratio")
parser$add_argument("-n", "--n_cells", type = "integer", help = "Number of cells (If specified, will override train_split)")
parser$add_argument("-m", "--method", choices = c("rff", "rff_gaussian", "pca", "scimilarity"), help = "Approximation method (None for full data/baseline)")
parser$add_argument("-t", "--task", choices = c("tissue", "celltype"), help = "")
parser$add_argument("-a", "--algorithm", default = "glmnet", choices = c("glmnet", "svm"), help = "Classification algorithm")
parser$add_argument("-D", "--D", type = "integer", help = "Number of RFF dimensions")
parser$add_argument("-C", "--n_components", type = "integer", help = "Number of PCA components")

args <- parser$parse_args()

if (!is.null(args$method)) {
  method <- args$method 
} else {
  method <- "baseline"      
}

replicate_id <- args$replicate

project_root <- "/scratch/zeynepaydin21/proj1"
data_dir <- file.path(project_root, "data", args$dataset)

source(file.path(project_root, "utils.R"))

get_peak_ram_gb <- function() {
  g <- gc()
  return(sum(g[, ncol(g)]) / 1024) 
}

if (args$dataset == "scea") {
  dataset_suffix <- paste0(args$dataset, "_", args$feature_mode)
} else {
  dataset_suffix <- args$dataset  # Other datasets always use intersection
}

# output directory structure:
# {project_root}/experiments/{dataset}/classification/{method}/{train_config}/rep{X}/

if (is.null(args$n_cells)) {
  output_dir <- file.path(project_root, "experiments", dataset_suffix, "classification", method, 
                          paste0("train", args$train_split*100), paste0("rep", replicate_id))
} else {
  output_dir <- file.path(project_root, "experiments", dataset_suffix, "classification", method, 
                          paste0("cells", args$n_cells), paste0("rep", replicate_id))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

########## PART 1: Preprocessing ##########
cat("Preprocessing data...\n")
preprocess_start <- Sys.time()
gc(reset = TRUE)

train_test_split <- function(data, y, train_split = args$train_split, n_cells = args$n_cells, seed = 2969 + replicate_id) {
  set.seed(seed)
  
  shuffle_idx <- sample(1:length(y))
  
  final_train_size <- if (!is.null(n_cells)) n_cells else floor(train_split * length(y))
  final_test_size <- if (!is.null(n_cells)) 5000 else length(y) - floor(train_split * length(y))
  
  if (final_test_size + final_train_size > length(y)) {
    stop(sprintf("Requested total size (train: %d + test: %d = %d) exceeds available cells (%d).", 
                 final_train_size, final_test_size, final_test_size + final_train_size, length(y)))
  }
  
  test_idx <- shuffle_idx[1:final_test_size]
  train_idx <- shuffle_idx[(final_test_size + 1):(final_test_size + final_train_size)]
  
  return(list(
    input_cells = nrow(data),
    input_genes = ncol(data),
    train_indices = train_idx,
    test_indices = test_idx,
    split_seed = seed
  ))
}

if (args$dataset == "scea") { 
  mtx <- qs_read(file.path(data_dir, paste0("scea_8tissues_", args$feature_mode, "_expression_matrix.qs")))
  labels <- qs_read(file.path(data_dir, paste0("scea_8tissues_", args$feature_mode, "_", args$task, "_labels.qs")))
} else {
  # For all other curated datasets, use intersection mode
  mtx <- qs_read(file.path(data_dir, paste0(args$dataset, "_intersection_expression_matrix.qs")))
  labels <- qs_read(file.path(data_dir, paste0(args$dataset, "_intersection_", args$task, "_labels.qs")))
}

cat("Creating train/test split...\n")
split_info <- train_test_split(data = mtx, y = labels)

X_train <- mtx[split_info$train_indices, , drop=FALSE]
X_test <- mtx[split_info$test_indices, , drop=FALSE]

all_labels_factor <- as.factor(labels)
y_train <- all_labels_factor[split_info$train_indices]
y_test <- all_labels_factor[split_info$test_indices]


rm(mtx, labels)

if (!(method %in% c("scimilarity"))) {
  cat("Normalizing and scaling datasets...\n")
  X_train <- preprocessor(X_train, scale_factor = 10000)
  X_test <- preprocessor(X_test, scale_factor = 10000)
}

preprocess_time <- as.numeric(difftime(Sys.time(), preprocess_start, units = "secs"))
preprocess_mem_gb <- get_peak_ram_gb()

########## PART 2: Approximation ##########
reduction_time <- 0
reduction_seed <- NULL
rff_metadata <- NULL
filename_base <- "baseline"

if (method != "baseline") {
  reduction_start <- Sys.time()
  
  if (method %in% c("rff", "rff_gaussian")) {

    if (method == "rff") {
      cat("Generating Random Fourier Features (Laplacian/Cauchy)...\n")
      kernel_type <- "laplacian"
      dist_method <- "manhattan"
      helper_script <- "calculate_sigma.R"

      generate_W <- function(n_features, n_components, param) {
        matrix(stats::rcauchy(n = n_features * n_components, location = 0, scale = param), 
               nrow = n_features, ncol = n_components)
      }
      
    } else { # rff_gaussian
      cat("Generating Random Fourier Features (Gaussian/Normal)...\n")
      kernel_type <- "gaussian"
      dist_method <- "euclidean"
      helper_script <- "calculate_sigma_gaussian.R"
      
      generate_W <- function(n_features, n_components, param) {
        matrix(stats::rnorm(n = n_features * n_components, mean = 0, sd = param), 
               nrow = n_features, ncol = n_components)
      }
    }
    
    D <- args$D
    
    if (nrow(X_train) >= 250000) {
      sigma_N <- 2500
    } else if (nrow(X_train) >= 100000) {
      sigma_N <- 2000
    } else {
      sigma_N <- 1000
    }
    
    reduction_seed <- 8876 + replicate_id
    set.seed(reduction_seed)
    
    distance_indices <- sample(1:nrow(X_train), sigma_N)
    sliced_data <- X_train[distance_indices, ]
    
    temp_dir <- tempdir()
    unique_id <- paste0(replicate_id, "_", Sys.getpid(), "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    temp_data_file <- file.path(temp_dir, paste0("data_", unique_id, ".qs"))
    temp_result_file <- file.path(temp_dir, paste0("sigma_", unique_id, ".qs"))
    qs_save(sliced_data, temp_data_file)
    
    cmd <- sprintf('Rscript --vanilla %s "%s" "%s"', 
                   helper_script, temp_data_file, temp_result_file)
    
    system_result <- system(cmd, wait = TRUE)
    
    sigma_results <- qs_read(temp_result_file)
    
    unlink(c(temp_data_file, temp_result_file))

    sigma <- sigma_results$sigma
    sigma_p <- sigma_results$sigma_p
    
    W <- generate_W(ncol(X_train), D, sigma_p)
    
    compute_rff <- function(X, W, D) {
      XW <- as.csr.matrix(X) %*% W
      Z <- sqrt(1 / D) * cbind(cos(XW), sin(XW))
      return(Z)
    }
    
    X_train <- compute_rff(X_train, W, D = D)
    X_test <- compute_rff(X_test, W, D = D)
    
    rff_metadata <- list(
      kernel = kernel_type,
      distance_method = dist_method,
      distance_indices = distance_indices,
      sigma = sigma,
      sigma_p = sigma_p,
      sigma_N = sigma_N
    )
    
    filename_base <- paste0("D", D, "_sigmaN", sigma_N)
    expected_dims <- 2 * D
    
  } else if (method == "pca") {
    cat("Performing PCA...\n")
    
    n_components <- args$n_components
    
    reduction_seed <- 752 + replicate_id
    set.seed(reduction_seed)
    
    svd_result <- RSpectra::svds(X_train, n_components, opts = list(center = TRUE, scale = FALSE))
    X_train <- svd_result$u %*% diag(svd_result$d)
    X_test <- X_test %*% svd_result$v 
    
    filename_base <- paste0("Ncomp", n_components)
    expected_dims <- n_components
    
  } else if (method == "scimilarity") {
    cat("Generating SCimilarity embeddings...\n")
    
    reduction_seed <- 4853 + replicate_id
    set.seed(reduction_seed)
    
    temp_dir <- tempdir()
    unique_id <- paste0(replicate_id, "_", Sys.getpid(), "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    indices_file <- file.path(temp_dir, paste0("scimilarity_indices_", unique_id, ".h5"))
    embeddings_file <- file.path(temp_dir, paste0("scimilarity_embeddings_", unique_id, ".h5"))

    # Save indices (convert to 0-based for Python)
    h5createFile(indices_file)
    h5write(as.integer(split_info$train_indices - 1), indices_file, "train_indices")
    h5write(as.integer(split_info$test_indices - 1), indices_file, "test_indices")
    H5close()
    
    python_script_path <- file.path(project_root, "get_scimilarity_embeddings.py")
    conda_path <- "module load conda3/latest && source /opt/ohpc/pub/compiler/conda3/latest/etc/profile.d/conda.sh && conda activate scimilarity"
      
    python_cmd <- sprintf(
      '%s && python "%s" --data_dir "%s" --feature_mode "%s" --indices_file "%s" --output_file "%s" --seed "%d"', 
      conda_path, python_script_path, data_dir, args$feature_mode, indices_file, embeddings_file, reduction_seed
    )
    
    python_output <- system(python_cmd, intern = TRUE)
    
    X_train <- t(h5read(embeddings_file, "train_embeddings"))
    X_test <- t(h5read(embeddings_file, "test_embeddings"))
    
    reduction_time <- h5readAttributes(embeddings_file, "/")$embedding_time
    preprocess_time <- h5readAttributes(embeddings_file, "/")$preprocess_time

    unlink(c(indices_file, embeddings_file))
    
    filename_base <- "scimilarity"
    expected_dims <- 128
  } 
  
  if (is.null(X_train) || is.null(X_test) || 
      ncol(X_train) != expected_dims || ncol(X_test) != expected_dims) {
    
    cat("ERROR: dimensionality reduction failed!\n")
    cat("Expected:", expected_dims, ", got Train:", ncol(X_train), "Test:", ncol(X_test), "\n")
    quit(status = 1)
  }
  
  
  cat("Dimensionality reduction successful:", nrow(X_train), "×", ncol(X_train), "\n")
  
  if (!(method %in% c("scimilarity")))  reduction_time <- as.numeric(difftime(Sys.time(), reduction_start, units = "secs"))
}

########## PART 3: Classification ##########

cat("Training", args$algorithm, "model...\n")

model_seed <- 65 + replicate_id
set.seed(model_seed)
modeling_mem_gb <- 0
gc(reset = TRUE)

if (args$algorithm == "glmnet") {
  model_start_time <- Sys.time()
  model <- glmnet::glmnet(x = X_train, 
                          y = y_train, 
                          family = "multinomial", 
                          alpha = 1)
  model_time <- as.numeric(difftime(Sys.time(), model_start_time, units = "secs"))
  
  best_lambda <- model$lambda[which.max(model$dev.ratio)]
  
  pred <- predict(model, newx = X_test, s = best_lambda, type = "class")
  pred <- factor(pred, levels = levels(y_test))
  
  prob <- predict(model, newx = X_test, s = best_lambda, type = "response")
  prob <- prob[, , 1]
  
} else if (args$algorithm == "svm") {
  model_start_time <- Sys.time()
  model <- e1071::svm(x = X_train, y = y_train, 
                      kernel = "linear", 
                      type = "C-classification", 
                      scale = TRUE, 
                      probability = TRUE)
  model_time <- as.numeric(difftime(Sys.time(), model_start_time, units = "secs"))
  
  pred <- predict(model, newdata = as.matrix(X_test))
  
  prob <- predict(model, newdata = as.matrix(X_test), probability = TRUE)
  prob <- attr(prob, "probabilities")
  
} 


conf_matrix <- caret::confusionMatrix(pred, y_test)
accuracy <- mean(pred == y_test)
precision <- mean(conf_matrix$byClass[, "Precision"], na.rm = TRUE)
recall <- mean(conf_matrix$byClass[, "Recall"], na.rm = TRUE)  
f1_score <- mean(conf_matrix$byClass[, "F1"], na.rm = TRUE)
auroc <- tryCatch({
  as.numeric(pROC::multiclass.roc(response = y_test, predictor = prob)$auc)
}, error = function(e) {
  NA
})

classification_script_time <- as.numeric(difftime(Sys.time(), preprocess_start, units = "secs"))

modeling_mem_gb <- get_peak_ram_gb()

train_sparsity <- as.list(calculate_sparsity(X_train))
test_sparsity <- as.list(calculate_sparsity(X_test))

metadata <- list(
  node_id = Sys.info()["nodename"],
  job_id = Sys.getenv("SLURM_JOB_ID"),
  args = args,
  method = method,
  split_info = split_info,
  scaling_method = "min_max",
  n_train = nrow(X_train),
  n_test = nrow(X_test),
  n_genes = ncol(X_test),
  n_train_classes = length(levels(y_train)),
  n_test_classes = length(levels(y_test)),
  train_sparsity = train_sparsity,
  test_sparsity = test_sparsity,
  
  reduction_seed = reduction_seed,
  model_seed = model_seed,
  rff_metadata = rff_metadata,
  
  metrics = list(
    accuracy = accuracy,
    precision = precision,
    recall = recall,
    f1_score = f1_score,
    auroc = auroc,
    confusion_matrix = conf_matrix
  ),
  
  predicted_classes = pred,
  class_probabilities = prob,
  # train_labels = y_train,
  true_labels = y_test,
  
  timing = list(
    preprocess_time = preprocess_time,
    reduction_time = reduction_time,
    model_training_time = model_time,
    classification_script_time = classification_script_time
  ),
  
  RAM_usage = list(
    preprocess_stage_peak_gb = round(preprocess_mem_gb, 3),
    modeling_R_peak_gb = round(modeling_mem_gb, 3)
  ),
  
  filename_base = filename_base
)

model_object <- list(
  trained_model = model,
  model_params = if (args$algorithm == "glmnet") {
    list(
      family = "multinomial",
      alpha = 1,
      best_lambda = best_lambda,
      lambda_selection = "max_dev_ratio"
    )
  } else if (args$algorithm == "svm") {
    list(
      kernel = "linear", 
      type = "C-classification",
      scale = TRUE,
      probability = TRUE
    )
  },
  training_context = list(
    r_version = R.version.string,
    package_versions = list(
      glmnet = if(args$algorithm == "glmnet") as.character(packageVersion("glmnet")) else NULL,
      e1071 = if(args$algorithm == "svm") as.character(packageVersion("e1071")) else NULL
    )
  )
)

metadata_filename <- paste0(filename_base, "_", args$algorithm, "_", args$task, "_metadata.RData")
model_filename <- paste0(filename_base, "_", args$algorithm, "_", args$task, "_model.RData")

save(metadata, file = file.path(output_dir, metadata_filename))
save(model_object, file = file.path(output_dir, model_filename))

cat("Script completed successfully!\n")