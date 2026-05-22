library(wordspace)

# wordspace::dist.matrix transposes M before passing to C++ code expecting dgCMatrix slots.
# With MatrixExtra loaded, t(dgCMatrix) -> dgRMatrix, breaking the C++ call.
# For sparse: route through TsparseMatrix (not overridden by MatrixExtra), then CsparseMatrix.
.prep_for_dist <- function(X) {
  if (inherits(X, "sparseMatrix")) {
    as(t(as(X, "TsparseMatrix")), "CsparseMatrix")
  } else {
    t(X)
  }
}

# Stratified sampling of sigma_N indices by class label.
# Equal allocation: floor(sigma_N / n_classes) per class, capped at class size.
stratified_sigma_indices <- function(labels, sigma_N) {
  labels <- as.factor(labels)
  classes <- levels(labels)
  sigma_N <- min(sigma_N, length(labels))

  class_idx <- lapply(classes, function(cl) which(labels == cl))
  class_sizes <- lengths(class_idx)

  n_per <- floor(sigma_N / length(classes))
  per_class <- pmin(n_per, class_sizes)

  unlist(mapply(function(idx, n) sample(idx, n),
                class_idx, per_class, SIMPLIFY = FALSE),
         use.names = FALSE)
}

estimate_sigma_laplacian <- function(X) {
  D <- wordspace::dist.matrix(.prep_for_dist(X), method = "manhattan", byrow = FALSE)
  sigma <- mean(D)
  sigma_p <- 1 / sigma
  list(sigma = sigma, sigma_p = sigma_p)
}

estimate_sigma_gaussian <- function(X) {
  D <- wordspace::dist.matrix(.prep_for_dist(X), method = "euclidean", byrow = FALSE)
  sigma <- mean(D)
  gamma <- 1 / (2 * sigma^2)
  sigma_p <- sqrt(2 * gamma)
  list(sigma = sigma, sigma_p = sigma_p)
}
