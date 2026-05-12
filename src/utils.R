library(Matrix)
library(sparseMatrixStats)

log_info <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

log_normalize <- function(X, scale_factor = 10000) {
  X <- as.csc.matrix(X)
  lib_sizes <- Matrix::rowSums(X)
  scaling_factors <- scale_factor / lib_sizes
  row_indices <- X@i + 1
  X@x <- X@x * scaling_factors[row_indices]
  log1p(X)
}

min_max_scale <- function(X, col_mins = NULL, col_maxs = NULL) {
  if (is.null(col_mins) || is.null(col_maxs)) {
    col_mins <- unname(sparseMatrixStats::colMins(X))
    col_maxs <- unname(sparseMatrixStats::colMaxs(X))
    if (any(col_mins > 0)) {
      warning(sprintf("SPARSITY ALERT: %d genes have non-zero minimums. Structural zeros in the test set will remain 0 for these genes.", sum(col_mins > 0)))
    }
  }
  denom <- col_maxs - col_mins
  denom[denom == 0] <- 1
  col_indices <- rep(seq_len(ncol(X)), diff(X@p))
  X@x <- (X@x - col_mins[col_indices]) / denom[col_indices]
  list(X = X, col_mins = col_mins, col_maxs = col_maxs)
}

select_hvg <- function(X_train, X_test, n_hvg) {
  gene_vars <- sparseMatrixStats::colVars(X_train)
  hvg_idx <- order(gene_vars, decreasing = TRUE)[1:n_hvg]
  list(
    X_train = X_train[, hvg_idx, drop = FALSE],
    X_test = X_test[, hvg_idx, drop = FALSE],
    hvg_indices = hvg_idx
  )
}

calculate_sparsity <- function(X) {
  total_elements <- as.numeric(nrow(X)) * as.numeric(ncol(X))
  non_zero_elements <- as.numeric(Matrix::nnzero(X))
  sparsity_ratio <- 1 - (non_zero_elements / total_elements)
  list(
    sparsity_ratio = sparsity_ratio,
    non_zero_elements = non_zero_elements,
    total_elements = total_elements
  )
}

get_peak_ram_gb <- function() {
  g <- gc()
  sum(g[, ncol(g)]) / 1024
}
