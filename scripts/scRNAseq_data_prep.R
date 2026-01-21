library(scRNAseq)
library(Matrix)
library(MatrixExtra)
library(qs2)
library(readr)
library(stringr)

create_data_dir <- function(dataset_name) {
  data_dir <- file.path(base_data_dir, dataset_name)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  return(data_dir)
}

filter_and_save <- function(mtx, labels, dataset_name, task, data_dir) {
  doublet_pattern <- "doublet|contaminant|contamination|low.quality|low_quality|unassigned|unknown|debris|multiplet"
  is_doublet <- grepl(doublet_pattern, labels, ignore.case = TRUE)
  
  cat("Removing", sum(is_doublet, na.rm = TRUE), "doublet/contaminant cells\n")
  
  valid_cells <- !is.na(labels) & !is_doublet
  mtx <- mtx[valid_cells, , drop=FALSE]
  labels <- as.factor(labels[valid_cells])
  
  cat("Matrix dimensions:", nrow(mtx), "cells x", ncol(mtx), "genes\n")
  cat("Number of classes:", length(levels(labels)), "\n")
  
  # saved matrix will be in dgr format since we transpose it every time
  matrix_file <- file.path(data_dir, paste0(dataset_name, "_intersection_expression_matrix.qs"))
  qs_save(mtx, matrix_file)
  
  label_file <- file.path(data_dir, paste0(dataset_name, "_intersection_", task, "_labels.qs"))
  qs_save(labels, label_file)
  cat("Save complete.")
}

ct_task <- "celltype"

# ensembl = TRUE not working properly


########## He 2020 ##########
gc()

# Human cortex single-nuclei RNA-seq data from He et al. (2017).
# Single-cell transcriptome profiling of an adult human cell atlas of 15 major organs. Genome Biol 21, 1:294.

sce <- HeOrganAtlasData() 
mtx <- t(assay(sce, "counts"))

tissue_labels <- sce$Tissue
ct_labels <- sce$Cell_type_in_merged_data

dataset <- "he_organs"
data_dir <- create_data_dir(dataset)

filter_and_save(mtx, ct_labels, dataset, ct_task, data_dir)
filter_and_save(mtx, tissue_labels, dataset, task = "tissue", data_dir)

########## Zhao 2020 ##########
rm(dataset, data_dir, sce, mtx, labels)
gc()

# Human liver immune single-cell RNA-seq data from Zhao et al. (2020).
# Single-cell RNA sequencing reveals the heterogeneity of liver-resident immune cells in human. Cell Discov 6, 22.

sce <- ZhaoImmuneLiverData()
mtx <- t(assay(sce, "counts"))
labels <- sce$fine
retained <- sce$retained

labels <- labels[retained]
mtx <- mtx[retained, , drop=FALSE]

dataset <- "zhao_immune"
data_dir <- create_data_dir(dataset)

filter_and_save(mtx, labels, dataset, ct_task, data_dir)

########## Zilionis 2019 ##########
rm(dataset, data_dir, sce, mtx, labels)
gc()

# Human/mouse lung cancer single-cell RNA-seq data from Zilionis et al. (2019).
# Single-cell transcriptomics of human and mouse lung cancers reveals conserved myeloid populations across individuals and species. Immunity 50(5), 1317-1334.

sce <- ZilionisLungData(which = "human") # default human
mtx <- t(assay(sce, "counts"))

labels <- sce$`Major cell type`
used <- sce$Used

labels <- labels[used]
mtx <- mtx[used, , drop=FALSE]

dataset <- "zilionis_lung"
data_dir <- create_data_dir(dataset)

filter_and_save(mtx, labels, dataset, ct_task, data_dir)