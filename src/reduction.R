.rff_gpu <- function(X_csr, W_mat, n_proj) {
  X_t <- torch:::torch__sparse_csr_tensor_unsafe(
    crow_indices = torch::torch_tensor(X_csr@p, dtype = torch::torch_int32()),
    col_indices = torch::torch_tensor(X_csr@j, dtype = torch::torch_int32()),
    values = torch::torch_tensor(X_csr@x, dtype = torch::torch_float32()),
    size = c(nrow(X_csr), ncol(X_csr))
  )$cuda()
  W_t <- torch::torch_tensor(W_mat, dtype = torch::torch_float32())$cuda()
  XW <- torch::torch_mm(X_t, W_t)
  Z <- sqrt(1 / n_proj) * torch::torch_cat(list(torch::torch_cos(XW), torch::torch_sin(XW)), dim = 2)
  as.matrix(Z$cpu())
}

rff_reduce <- function(X_train, X_test, n_dim, run_id, kernel = "laplacian",
                       use_gpu = FALSE) {
  t_start <- Sys.time()

  if (kernel == "laplacian") {
    log_info("Generating Random Fourier Features (Laplacian/Cauchy)...")
    dist_method <- "manhattan"
    generate_W <- function(n_features, n_components, param) {
      matrix(stats::rcauchy(n = n_features * n_components, location = 0, scale = param),
             nrow = n_features, ncol = n_components)
    }
  } else {
    log_info("Generating Random Fourier Features (Gaussian/Normal)...")
    dist_method <- "euclidean"
    generate_W <- function(n_features, n_components, param) {
      matrix(stats::rnorm(n = n_features * n_components, mean = 0, sd = param),
             nrow = n_features, ncol = n_components)
    }
  }

  if (n_dim %% 2 != 0) stop("--n_dim must be even for RFF (output = 2 * n_projections)")
  n_proj <- n_dim / 2L

  sigma_N <- 500
  reduction_seed <- 8876 + run_id
  set.seed(reduction_seed)

  distance_indices <- sample(1:nrow(X_train), sigma_N)
  sigma_results <- if (kernel == "laplacian") {
    estimate_sigma_laplacian(X_train[distance_indices, ])
  } else {
    estimate_sigma_gaussian(X_train[distance_indices, ])
  }
  sigma <- sigma_results$sigma
  sigma_p <- sigma_results$sigma_p

  W <- generate_W(ncol(X_train), n_proj, sigma_p)

  compute_rff <- function(X, W, n_proj) {
    X_csr <- as.csr.matrix(X)
    if (use_gpu) {
      .rff_gpu(X_csr, W, n_proj)
    } else {
      XW <- X_csr %*% W
      sqrt(1 / n_proj) * cbind(cos(XW), sin(XW))
    }
  }

  list(
    X_train = compute_rff(X_train, W, n_proj),
    X_test = compute_rff(X_test, W, n_proj),
    rff_metadata = list(
      kernel = kernel,
      distance_method = dist_method,
      distance_indices = distance_indices,
      sigma = sigma,
      sigma_p = sigma_p,
      sigma_N = sigma_N
    ),
    filename_base = paste0("dim", n_dim, "_sigmaN", sigma_N),
    expected_dims = n_dim,
    reduction_seed = reduction_seed,
    reduction_time = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  )
}

pca_reduce <- function(X_train, X_test, n_dim, run_id) {
  log_info("Performing PCA...")
  t_start <- Sys.time()

  reduction_seed <- 752 + run_id
  set.seed(reduction_seed)

  train_means <- sparseMatrixStats::colMeans2(X_train)
  svd_result <- RSpectra::svds(X_train, n_dim, opts = list(center = TRUE, scale = FALSE))
  X_train_r <- svd_result$u %*% diag(svd_result$d)
  mean_shift <- matrix(train_means, nrow = 1) %*% svd_result$v
  X_test_r <- sweep(X_test %*% svd_result$v, 2, as.vector(mean_shift), FUN = "-")

  list(
    X_train = X_train_r,
    X_test = X_test_r,
    rff_metadata = NULL,
    filename_base = paste0("dim", n_dim),
    expected_dims = n_dim,
    reduction_seed = reduction_seed,
    reduction_time = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  )
}

scvi_embed <- function(split_info, data_dir, project_root, run_id) {
  log_info("Generating scVI (pretrained Census) embeddings...")
  suppressPackageStartupMessages(library(rhdf5))

  temp_dir <- tempdir()
  unique_id <- paste0(run_id, "_", Sys.getpid(), "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  indices_file    <- file.path(temp_dir, paste0("scvi_indices_", unique_id, ".h5"))
  embeddings_file <- file.path(temp_dir, paste0("scvi_embeddings_", unique_id, ".h5"))

  h5createFile(indices_file)
  h5write(as.integer(split_info$train_indices - 1), indices_file, "train_indices")
  h5write(as.integer(split_info$test_indices - 1), indices_file, "test_indices")
  H5close()

  python_script_path <- file.path(project_root, "src", "get_scvi_embeddings.py")
  gene_map_file <- file.path(project_root, "data", "gene_sets", "ensembl_symbol_to_ensembl_id.tsv")
  python_cmd <- sprintf(
    '$HOME/.conda/envs/scvi/bin/python "%s" --data_dir "%s" --indices_file "%s" --output_file "%s" --gene_map_file "%s" 2>&1',
    python_script_path, data_dir, indices_file, embeddings_file, gene_map_file
  )
  exit_code <- system(python_cmd)
  if (exit_code != 0) stop("scVI Python script failed (exit code ", exit_code, ")")

  X_train_r <- t(h5read(embeddings_file, "train_embeddings"))
  X_test_r <- t(h5read(embeddings_file, "test_embeddings"))
  attrs <- h5readAttributes(embeddings_file, "/")
  unlink(c(indices_file, embeddings_file))

  list(
    X_train = X_train_r,
    X_test = X_test_r,
    rff_metadata = NULL,
    filename_base = "scvi",
    expected_dims = 50L,
    reduction_seed = NA_integer_,
    reduction_time = attrs$embedding_time,
    preprocess_time = attrs$preprocess_time
  )
}

scimilarity_embed <- function(split_info, data_dir, project_root, run_id) {
  log_info("Generating SCimilarity embeddings...")
  suppressPackageStartupMessages(library(rhdf5))

  reduction_seed <- 4853 + run_id
  temp_dir <- tempdir()
  unique_id <- paste0(run_id, "_", Sys.getpid(), "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  indices_file    <- file.path(temp_dir, paste0("scimilarity_indices_", unique_id, ".h5"))
  embeddings_file <- file.path(temp_dir, paste0("scimilarity_embeddings_", unique_id, ".h5"))

  h5createFile(indices_file)
  h5write(as.integer(split_info$train_indices - 1), indices_file, "train_indices")
  h5write(as.integer(split_info$test_indices - 1), indices_file, "test_indices")
  H5close()

  python_script_path <- file.path(project_root, "src", "get_scimilarity_embeddings.py")
  python_cmd <- sprintf(
    '$HOME/.conda/envs/scimilarity/bin/python "%s" --data_dir "%s" --indices_file "%s" --output_file "%s" --seed "%d" 2>&1',
    python_script_path, data_dir, indices_file, embeddings_file, reduction_seed
  )
  exit_code <- system(python_cmd)
  if (exit_code != 0) stop("SCimilarity Python script failed (exit code ", exit_code, ")")

  X_train_r <- t(h5read(embeddings_file, "train_embeddings"))
  X_test_r <- t(h5read(embeddings_file, "test_embeddings"))
  attrs <- h5readAttributes(embeddings_file, "/")
  unlink(c(indices_file, embeddings_file))

  list(
    X_train = X_train_r,
    X_test = X_test_r,
    rff_metadata = NULL,
    filename_base = "scimilarity",
    expected_dims = 128,
    reduction_seed = reduction_seed,
    reduction_time = attrs$embedding_time,
    preprocess_time = attrs$preprocess_time
  )
}
