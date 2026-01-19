suppressPackageStartupMessages({
  library(Matrix)
  library(wordspace)
  library(qs2)
})

gc(reset = TRUE)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

sliced_data <- qs_read(input_file)

pdist <- function(X1, X2, method = "manhattan") {
  if (identical(X1, X2) == TRUE) {
    D <- wordspace::dist.matrix(X1, method = method)
  } else {
    D <- wordspace::dist.matrix(rbind(X1, X2), method = method)
    D <- D[1:nrow(X1), (nrow(X1) + 1):(nrow(X1) + nrow(X2))]
  }
  return(D)
}

D_train <- pdist(sliced_data, sliced_data)

sigma <- mean(D_train)
gamma <- 1 / sigma
sigma_p <- gamma

final_gc <- gc()
mb_val <- which(colnames(final_gc) == "max used") + 1 # MB
peak_memory_gb <- max(final_gc[, mb_val]) / 1024

result <- list(
  sigma = sigma,
  sigma_p = sigma_p,
  peak_memory_gb = peak_memory_gb
)

qs_save(result, output_file)