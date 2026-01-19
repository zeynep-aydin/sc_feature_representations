library(Matrix)
library(sparseMatrixStats)

#' Preprocess sparse matrix with log normalization and min-max scaling
#'
#' @param X Sparse matrix (cells x genes)
#' @param scale_factor Scaling factor for library size normalization (default: 10000)
#' @return Preprocessed Csparse matrix
preprocessor <- function(X, scale_factor = 10000) {
  X <- as.csc.matrix(X) # back to CSC format

  # log normalizing
  lib_sizes <- Matrix::rowSums(X)
  scaling_factors <- scale_factor / lib_sizes

  row_indices <- X@i + 1
  X@x <- X@x * scaling_factors[row_indices]
  X <- log1p(X)

  # scaling
  col_mins <- unname(sparseMatrixStats::colMins(X))
  col_maxs <- unname(sparseMatrixStats::colMaxs(X))
  denom <- col_maxs - col_mins
  denom[denom == 0] <- 1

  col_indices <- rep(seq_len(ncol(X)), diff(X@p))
  X@x <- (X@x - col_mins[col_indices]) / denom[col_indices]
  return(X)
}


#' Calculate sparsity of a sparse matrix
#'
#' @param X Sparse matrix
#' @return List with sparsity statistics
calculate_sparsity <- function(X) {
  total_elements <- as.numeric(nrow(X)) * as.numeric(ncol(X))
  non_zero_elements <- as.numeric(Matrix::nnzero(X))
  sparsity_ratio <- 1 - (non_zero_elements / total_elements)
  return(list(
    sparsity_ratio = sparsity_ratio,
    non_zero_elements = non_zero_elements,
    total_elements = total_elements
  ))
}