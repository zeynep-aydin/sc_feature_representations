#!/usr/bin/env Rscript
#
# Prepare CRC anti-PD-1 (GSE236581, Cancer Cell 2024) for the benchmark pipeline.
# Tasks: cell type (MajorCellType), tissue state (Tumor/Normal), compound (tissue_celltype).
#
# Source: GEO supplementary files -- gene-symbol counts (genes x cells) + per-cell metadata.
#   counts: GSE236581_counts.mtx.gz, GSE236581_barcodes.tsv.gz, GSE236581_features.tsv.gz
#   metadata: GSE236581_CRC-ICB_metadata.txt.gz
# Subset to Tissue %in% c("Tumor", "Normal") -- Blood/LN/TN excluded, not part of the
# tumor-vs-normal / cell-type classification tasks.

suppressPackageStartupMessages({
  library(Matrix)
  library(MatrixExtra)
})

script_dir <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "utils_prep.R"))

PROJ_ROOT <- Sys.getenv("PROJ_ROOT", unset = dirname(script_dir))
raw_dir <- file.path(PROJ_ROOT, "raw_data", "crc_icb")
out_dir <- file.path(PROJ_ROOT, "data", "crc_icb")

cat("Reading genes/barcodes...\n")
genes <- read.table(gzfile(file.path(raw_dir, "GSE236581_features.tsv.gz")),
                     header = FALSE, sep = "\t", stringsAsFactors = FALSE)[[1]]  # col1 = col2 = symbol
barcodes <- readLines(gzcon(file(file.path(raw_dir, "GSE236581_barcodes.tsv.gz"), "rb")))
cat(sprintf("  %d genes | %d barcodes\n", length(genes), length(barcodes)))

cat("Reading matrix.mtx (~4GB, may take a few minutes)...\n")
mtx <- readMM(gzfile(file.path(raw_dir, "GSE236581_counts.mtx.gz")))
mtx <- t(mtx)
mtx <- as(mtx, "CsparseMatrix")
rownames(mtx) <- barcodes
colnames(mtx) <- genes
cat(sprintf("  matrix: %d cells x %d genes\n", nrow(mtx), ncol(mtx)))

cat("Reading metadata...\n")
meta <- read.table(gzfile(file.path(raw_dir, "GSE236581_CRC-ICB_metadata.txt.gz")),
                    header = TRUE, sep = " ", quote = "\"", row.names = 1, check.names = FALSE)

common <- intersect(rownames(mtx), rownames(meta))
cat(sprintf("Cells with metadata: %d / %d\n", length(common), nrow(mtx)))
mtx <- mtx[common, , drop = FALSE]
meta <- meta[common, , drop = FALSE]

ct <- meta$MajorCellType
tissue <- meta$Tissue
patient <- meta$Patient

# Subset to Tumor + Normal (exclude Blood/LN/TN -- not part of the classification tasks)
keep <- tissue %in% c("Tumor", "Normal")
cat(sprintf("Keeping %d / %d Tumor+Normal cells (dropping %s)\n",
            sum(keep), length(keep), paste(sort(unique(tissue[!keep])), collapse = ", ")))
mtx <- mtx[keep, , drop = FALSE]
ct <- ct[keep]; tissue <- tissue[keep]; patient <- patient[keep]

compound <- paste0(tolower(tissue), "_", ct)

d <- drop_by_min_cells(mtx, ct, patient, min_cells = 50, extra = list(tissue = tissue, compound = compound))
mtx <- d$mtx; ct <- d$ct; patient <- d$don; tissue <- d$tissue; compound <- d$compound

cat("\nCell type (MajorCellType) counts:\n"); print(sort(table(ct), decreasing = TRUE))
cat("\nTissue counts:\n"); print(table(tissue))
cat("\nCompound (tissue_celltype) counts:\n"); print(sort(table(compound), decreasing = TRUE))
cat("\nDonor counts:\n"); print(sort(table(patient), decreasing = TRUE))

save_dataset(mtx, ct, patient, out_dir)
qs2::qs_save(as.factor(tissue), file.path(out_dir, "tissue.qs"))
qs2::qs_save(droplevels(as.factor(compound)), file.path(out_dir, "labels_compound.qs"))
export_h5ad(mtx, out_dir, script_dir)

cat("\nDone.\n")
