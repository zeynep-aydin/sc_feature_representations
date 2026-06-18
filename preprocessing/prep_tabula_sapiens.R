#!/usr/bin/env Rscript
#
# Prepare Tabula Sapiens for the benchmark pipeline.
#
# Reads filtered data/tabula_sapiens/counts.h5ad (written by prep_tabula_sapiens.py)
# and saves counts.qs / labels.qs / donors.qs / tissue.qs.

suppressPackageStartupMessages({
  library(Matrix)
  library(rhdf5)
})

script_dir <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "utils_prep.R"))

PROJ_ROOT <- Sys.getenv("PROJ_ROOT", unset = dirname(script_dir))
out_dir <- file.path(PROJ_ROOT, "data", "tabula_sapiens")
h5ad <- file.path(out_dir, "counts.h5ad")

read_h5ad_col <- function(path, col) {
  node <- paste0("/obs/", col)
  tryCatch({
    cats  <- as.character(h5read(path, paste0(node, "/categories")))
    codes <- as.integer(h5read(path, paste0(node, "/codes"))) + 1L
    cats[codes]
  }, error = function(e) {
    as.character(h5read(path, node))
  })
}

cat("Reading genes and obs...\n")
genes <- as.character(h5read(h5ad, "/var/_index"))
cell_types <- read_h5ad_col(h5ad, "broad_cell_class")
cell_types_fine <- read_h5ad_col(h5ad, "cell_type")
compartments <- read_h5ad_col(h5ad, "compartment")
donors <- read_h5ad_col(h5ad, "donor_id")
tissues <- read_h5ad_col(h5ad, "tissue_in_publication")
n_cells <- length(cell_types)
n_genes <- length(genes)
cat(sprintf("  %d cells x %d genes\n", n_cells, n_genes))

cat("Reading sparse matrix...\n")
x <- as.numeric(h5read(h5ad, "X/data"))
indices <- as.integer(h5read(h5ad, "X/indices"))
indptr <- as.integer(h5read(h5ad, "X/indptr"))

# h5ad stores CSR (cells x genes); dgCMatrix is CSC, build as (genes x cells) then transpose
mtx <- t(new("dgCMatrix",
             x = x,
             i = indices,
             p = indptr,
             Dim = c(n_genes, n_cells),
             Dimnames = list(genes, NULL)))
validObject(mtx)
rm(x, indices, indptr); gc()
cat(sprintf("Assembled: %d x %d, NNZ=%d\n", nrow(mtx), ncol(mtx), nnzero(mtx)))

save_dataset(mtx, cell_types, donors, out_dir)
qs2::qs_save(droplevels(as.factor(cell_types_fine)), file.path(out_dir, "labels_fine.qs"))
qs2::qs_save(droplevels(as.factor(compartments)), file.path(out_dir, "labels_compartment.qs"))
qs2::qs_save(as.factor(tissues), file.path(out_dir, "tissue.qs"))

cat("\nBroad cell class counts:\n"); print(table(cell_types))
cat("\nCompartment counts:\n"); print(table(compartments))
cat("\nTissue counts:\n"); print(table(tissues))
