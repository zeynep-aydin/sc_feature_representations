#!/usr/bin/env Rscript
#
# Prepare zilionis_lung for the benchmark pipeline.
#
# Source: ZilionisLungData(which = "human") from the scRNAseq Bioconductor package.
# Filters: authors' Used flag, then label-pattern drop (doublets + Patient*-specific + ND).
# Output:  data/zilionis_lung/{counts,labels,donors}.qs + counts.h5ad

suppressPackageStartupMessages({
  library(scRNAseq)
  library(SingleCellExperiment)
  library(MatrixExtra)
})

script_dir <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "utils_prep.R"))

PROJ_ROOT <- Sys.getenv("PROJ_ROOT", unset = dirname(script_dir))
out_dir   <- file.path(PROJ_ROOT, "data", "zilionis_lung")

cat("Loading ZilionisLungData\n")
sce <- ZilionisLungData(which = "human")

mtx <- t(assay(sce, "counts"))  # cells x genes
ct  <- sce$`Major cell type`
don <- sce$Patient
tissue <- sce$Tissue
used <- sce$Used

cat(sprintf("Raw: %d cells, %d genes\n", nrow(mtx), ncol(mtx)))

mtx <- mtx[used, , drop = FALSE]
ct  <- ct[used]
don <- don[used]
tissue <- tissue[used]
cat(sprintf("After Used filter: %d cells\n", nrow(mtx)))

# Drop doublets + Patient*-specific (donor-exclusive types) + ND (undefined)
d <- drop_by_label(mtx, ct, don, extra_pattern = "Patient.*specific|^ND$",
                   extra = list(tissue = tissue))
mtx <- d$mtx; ct <- d$ct; don <- d$don; tissue <- d$tissue

d <- drop_by_min_cells(mtx, ct, don, min_cells = 50, extra = list(tissue = tissue))
mtx <- d$mtx; ct <- d$ct; don <- d$don; tissue <- d$tissue

cat("\nFinal cell type counts:\n")
print(sort(table(ct), decreasing = TRUE))
cat("\nDonor counts:\n")
print(sort(table(don), decreasing = TRUE))
cat("\nTissue counts:\n")
print(table(tissue))

save_dataset(mtx, ct, don, out_dir)
qs2::qs_save(as.factor(tissue), file.path(out_dir, "tissue.qs"))
export_h5ad(mtx, out_dir, script_dir)
