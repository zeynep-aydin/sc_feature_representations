#!/usr/bin/env Rscript
#
# Prepare Lodi et al. 2025 pan-cancer TME atlas (Cell Reports Medicine) for the
# benchmark pipeline. Task: 9-way cancer-type-of-origin classification.
#
# Source: lab portal MTX triplet (raw integer counts, gene symbols) + per-cell
# metadata CSV. Downloaded to raw_data/pancancer/.
#   counts: Lodi2025_counts/{matrix.mtx,genes.tsv,barcodes.tsv}  (genes x cells)
#   metadata: 9209-Lodi2025_metadata.csv.zip -> Lodi2025_metadata_update.csv
#
# Produces two datasets:
#   data/lodi_pancancer/     malignant only (Cancer/epithelial, Normal site excluded) -- primary
#   data/lodi_pancancer_all/ all cell types (Normal site excluded)                    -- ablation
# In both: labels = TumorType (cancer type), donors = Patient.

suppressPackageStartupMessages({
  library(Matrix)
  library(MatrixExtra)
})

script_dir <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "utils_prep.R"))

PROJ_ROOT <- Sys.getenv("PROJ_ROOT", unset = dirname(script_dir))
raw_dir <- file.path(PROJ_ROOT, "raw_data", "lodi_pancancer")
counts_dir <- file.path(raw_dir, "Lodi2025_counts")
meta_zip <- file.path(raw_dir, "9209-Lodi2025_metadata.csv.zip")

# Read counts
cat("Reading genes/barcodes...\n")
genes    <- read.table(file.path(counts_dir, "genes.tsv"), header = FALSE,
                       sep = "\t", stringsAsFactors = FALSE)[[1]]  # col 1 = gene symbol (unique)
barcodes <- readLines(file.path(counts_dir, "barcodes.tsv"))
cat(sprintf("  %d genes | %d barcodes\n", length(genes), length(barcodes)))

cat("Reading matrix.mtx\n")
mtx <- readMM(file.path(counts_dir, "matrix.mtx"))
mtx <- t(mtx)                                        
mtx <- as(mtx, "CsparseMatrix")
rownames(mtx) <- barcodes
colnames(mtx) <- genes
cat(sprintf("  matrix: %d cells x %d genes\n", nrow(mtx), ncol(mtx)))

# Read metadata and align to matrix rows
cat("Reading metadata...\n")
meta <- read.csv(unz(meta_zip, "Lodi2025_metadata_update.csv"), stringsAsFactors = FALSE)
rownames(meta) <- meta$Barcodes

common <- intersect(rownames(mtx), rownames(meta))
cat(sprintf("Cells with metadata: %d / %d\n", length(common), nrow(mtx)))
mtx  <- mtx[common, , drop = FALSE]
meta <- meta[common, , drop = FALSE]

tumortype <- meta$TumorType
patient   <- meta$Patient
majorct   <- meta$Majorcelltype_annotation
biopsy    <- meta$BiopsySite

# Drop cancer types that can't be donor-split (<2 patients)
# 'endometrial' has a single patient -> drop. 
# (Donor-based train/test needs >=1 donor per class on each side.)
don_per_type <- tapply(patient, tumortype, function(x) length(unique(x)))
drop_types <- names(don_per_type[don_per_type < 2])
if (length(drop_types) > 0) {
  cat(sprintf("Dropping %d type(s) with <2 patients: %s\n",
              length(drop_types), paste(drop_types, collapse = ", ")))
  keep <- !tumortype %in% drop_types
  mtx <- mtx[keep, , drop = FALSE]
  tumortype <- tumortype[keep]; patient <- patient[keep]
  majorct <- majorct[keep]; biopsy <- biopsy[keep]
}

# Exclude Normal-biopsy cells
# Cancer/epithelial cells from Normal sites are normal epithelium, mislabeled
# into a cancer class. Keep Tumor + LNt (lymph-node tumor) + Metastasis.
keep <- biopsy != "Normal"
cat(sprintf("Dropping %d Normal-biopsy cells (keeping %s)\n",
            sum(!keep), paste(sort(unique(biopsy[keep])), collapse = ", ")))
mtx <- mtx[keep, , drop = FALSE]
tumortype <- tumortype[keep]; patient <- patient[keep]
majorct <- majorct[keep]; biopsy <- biopsy[keep]

save_variant <- function(sub_mtx, sub_tt, sub_don, out_name) {
  out_dir <- file.path(PROJ_ROOT, "data", out_name)
  d <- drop_by_min_cells(sub_mtx, sub_tt, sub_don, min_cells = 50)
  cat(sprintf("\n[%s] cancer-type counts:\n", out_name))
  print(sort(table(d$ct), decreasing = TRUE))
  cat("[", out_name, "] patients per type:\n", sep = "")
  print(tapply(d$don, d$ct, function(x) length(unique(x))))
  save_dataset(d$mtx, d$ct, d$don, out_dir)
  export_h5ad(d$mtx, out_dir, script_dir)
}

# Primary: malignant only
mal <- majorct == "Cancer/epithelial"
cat(sprintf("\nMalignant (Cancer/epithelial) cells: %d / %d\n", sum(mal), length(mal)))
save_variant(mtx[mal, , drop = FALSE], tumortype[mal], patient[mal], "lodi_pancancer")

# Ablation: all cell types
save_variant(mtx, tumortype, patient, "lodi_pancancer_all")

cat("\nDone.\n")
