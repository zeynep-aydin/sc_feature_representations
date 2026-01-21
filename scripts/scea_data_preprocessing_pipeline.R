# Directory structure
# ├── cache/                            # Reusable across feature modes
# │   ├── participants.rds            
# │   ├── cell_mapping_filtered.rds   
# │   ├── common_features_thresh08.rds
# │   ├── common_features_intersection.rds
# │   ├── cell_mapping_E-ANND-1.rds  	# Individual cohort mappings
# │   └── ...
# ├── cohort_specific/                  # Feature-mode specific matrices
# │   ├── MTX_E-ANND-1_feature_mode_threshold.qs
# │   └── ...
# └── cleaned_data/                     # Final datasets
#     ├── logs/
#     ├── scea_8tissues_thresh08_expression_matrix.qs
#     ├── scea_8tissues_thresh08_tissue_labels.qs 
#     ├── scea_8tissues_thresh08_celltype_labels.qs
#     ├── scea_8tissues_thresh08_barcodes.qs
#     ├── scea_8tissues_thresh08_features.qs
#     ├── scea_8tissues_thresh08_metadata.rds
#     ├── scea_8tissues_intersection_expression_matrix.qs
#     └── scea_8tissues_intersection_tissue_labels.qs
#     └── ...

########## Config ########## 
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(Matrix)
  library(MatrixExtra)
  library(qs2)
})


CONFIG <- list(
  assay_dir = file.path("raw_data", "assays"),
  metadata_dir = file.path("raw_data", "metadatas"),
  output_dir = "cleaned_data",
  cache_dir = file.path("cleaned_data", "cache"),
  cohort_specific = file.path("cleaned_data", "cohort_specific"),
  tissues = list(
    kidney = "kidney",
    liver = "liver", 
    lung = "lung",
    colon = "colon",
    heart = "heart",
    pancreas = "pancreas",
    brain = "brain",
    prostate = "prostate"
  ),
  min_cells = 200, # minimum cells per cell type per tissue
  feature_mode = "intersection", # either "intersection" or "threshold"
  feature_threshold = 0.8, # used only in threshold mode, e.g. features must be present in 80% of cohorts per tissue
  tissue_exclusions = list(
    lung = c("lung bud", "lung-draining lymph node")
  ),
  force_rerun = list(
    participants = FALSE,       
    features = TRUE,            
    cell_mapping = FALSE,       
    matrices = TRUE             
  )
)


dir.create(CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(CONFIG$output_dir, "logs"), showWarnings = FALSE)
dir.create(CONFIG$cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(CONFIG$cohort_specific, recursive = TRUE, showWarnings = FALSE)

log_info <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] %s: %s", timestamp, level, message)
  cat(log_entry, "\n")

  write(log_entry, file.path(CONFIG$output_dir, "logs", "pipeline.log"), append = TRUE)
}

cache_exists <- function(cache_name) {
  cache_file <- file.path(CONFIG$cache_dir, paste0(cache_name, ".rds"))
  file.exists(cache_file)
}

save_cache <- function(object, cache_name) {
  cache_file <- file.path(CONFIG$cache_dir, paste0(cache_name, ".rds"))
  saveRDS(object, cache_file)
  log_info(paste("Cached result saved:", cache_name))
}

load_cache <- function(cache_name) {
  cache_file <- file.path(CONFIG$cache_dir, paste0(cache_name, ".rds"))
  log_info(paste("Loading cached result:", cache_name))
  readRDS(cache_file)
}

# Generate dataset name based on config
generate_dataset_name <- function() {
  tissue_count <- length(CONFIG$tissues)
  mode_str <- if (CONFIG$feature_mode == "threshold") {
    paste0("thresh", gsub("\\.", "", as.character(CONFIG$feature_threshold)))
  } else {
    CONFIG$feature_mode
  }
  paste0("scea_", tissue_count, "tissues_", mode_str)
}

########## Identifying participants with target tissues ########## 
identify_participants <- function() {
  cache_name <- "participants"
  
  # Check if cached result exists and we don't need to force rerun
  if (cache_exists(cache_name) && !CONFIG$force_rerun$participants) {
    log_info("=== STEP 1: Loading cached participants ===")
    return(load_cache(cache_name))
  }
  
  log_info("=== Step 1: Identifying participants with target tissues ===")
  
  sdrf_files <- list.files(CONFIG$metadata_dir, pattern = "\\.sdrf.txt$", full.names = TRUE)
  participants <- list()
  
  for (sdrf_file in sdrf_files) {
    cohort <- sub("__.*", "", basename(sdrf_file))
    log_info(paste("Processing cohort:", cohort))
    
    # Read sdrf file for sample metadata
    sdrf <- read_tsv(sdrf_file, col_types = cols(.default = "c"), show_col_types = FALSE)
    colnames(sdrf) <- tolower(gsub(" ", "_", colnames(sdrf)))
    
    # Read cell metadata file for cell type info
    cell_meta_file <- file.path(CONFIG$metadata_dir, paste0(cohort, "__", cohort, ".cell_metadata.tsv"))
    cmdf <- read_tsv(cell_meta_file, col_types = cols(.default = "c"), show_col_types = FALSE)
    
    # Find relevant columns
    author_cell_type_col <- grep("inferred.*cell.*type.*authors.*labels?", colnames(cmdf), value = TRUE)[1]
    ontology_cell_type_col <- grep("inferred.*cell.*type.*ontology.*labels?", colnames(cmdf), value = TRUE)[1]
    organism_col <- grep("\\borganism\\b", colnames(sdrf), value = TRUE)[1]
    individual_col <- grep("\\bindividual\\b", colnames(sdrf), value = TRUE)[1] 
    tissue_col <- grep("organism[\\. _-]*part", colnames(sdrf), value = TRUE)[1]
    if (any(is.na(c(organism_col, individual_col, tissue_col)))) {
      log_info(paste("Skipping", cohort, "--- missing required columns"))
      next
    } else if (all(is.na(c(author_cell_type_col, ontology_cell_type_col)))) {
      log_info(paste("Skipping", cohort, "--- missing cell type columns"))
      next
    } 
    
    # Filter for human samples
    human_data <- sdrf %>%
      filter(.data[[organism_col]] == "Homo sapiens") %>%
      filter(!is.na(.data[[individual_col]]), !is.na(.data[[tissue_col]]))
    
    if (nrow(human_data) == 0) next
    
    # Extract participants for each tissue
    participants[[cohort]] <- list()
    for (tissue in names(CONFIG$tissues)) {
      tissue_pattern <- CONFIG$tissues[[tissue]]
      
      # Apply tissue-specific exclusions if they exist
      exclusions <- CONFIG$tissue_exclusions[[tissue]]
      
      tissue_participants <- human_data %>%
        filter(str_detect(tolower(.data[[tissue_col]]), tissue_pattern))
      
      if (!is.null(exclusions) && length(exclusions) > 0) {
        tissue_participants <- tissue_participants %>%
          filter(!str_detect(tolower(.data[[tissue_col]]), paste(exclusions, collapse = "|")))
      }
      
      tissue_participants <- tissue_participants %>%
        pull(.data[[individual_col]]) %>%
        unique()
      
      if (length(tissue_participants) > 0) {
        participants[[cohort]][[tissue]] <- tissue_participants
      }
    }
    
    # Remove cohorts with no participants
    if (length(participants[[cohort]]) == 0) participants[[cohort]] <- NULL
  }
  
  save_cache(participants, cache_name)
  # saveRDS(participants, file.path(CONFIG$output_dir, "participants.rds"))
  log_info(paste("Found participants in", length(participants), "cohorts"))
  return(participants)
}

########## Find common features across cohorts ########## 
find_common_features <- function(participants) {
  cache_name <- paste0("common_features_", CONFIG$feature_mode)
  if (CONFIG$feature_mode == "threshold") {
    cache_name <- paste0(cache_name, "_", gsub("\\.", "", as.character(CONFIG$feature_threshold)))
  }
  
  if (cache_exists(cache_name) && !CONFIG$force_rerun$features) {
    log_info("=== STEP 2: Loading cached common features ===")
    result_obj <- load_cache(cache_name)
    result <- result_obj$global
    log_info(paste("Global intersection:", length(result), "features"))
    return(result)
  }
  
  log_info("=== Step 2: Finding common features across cohorts ===")
  
  cohorts <- names(participants)
  features_by_tissue <- setNames(vector("list", length(CONFIG$tissues)), names(CONFIG$tissues))
  
  for (tissue in names(CONFIG$tissues)) {
    features_by_tissue[[tissue]] <- list()
  }
  
  for (cohort in cohorts) {
    ftr_file <- file.path(CONFIG$assay_dir, cohort, 
                          paste0(cohort, ".aggregated_filtered_normalised_counts.mtx_rows.gz"))
    
    features <- read_tsv(ftr_file, col_names = "ensembl_id", show_col_types = FALSE)$ensembl_id
    
    # Add features for each tissue this cohort has
    for (tissue in names(participants[[cohort]])) {
      if (length(participants[[cohort]][[tissue]]) > 0) {
        features_by_tissue[[tissue]][[cohort]] <- features
      }
    }
  }
  
  # Find intersection of features for each tissue
  common_features_by_tissue <- list()
  
  for (tissue in names(features_by_tissue)) {
    tissue_cohorts <- features_by_tissue[[tissue]]
    
    if (length(tissue_cohorts) == 1) {
      common_features_by_tissue[[tissue]] <- tissue_cohorts[[1]]
    } else {
      if (CONFIG$feature_mode == "threshold") {
        # Thresholding mode: keep features present in >= threshold * cohorts
        all_features <- unlist(tissue_cohorts)
        feature_counts <- table(all_features)
        min_cohorts <- floor(CONFIG$feature_threshold * length(tissue_cohorts))
        common_features_by_tissue[[tissue]] <- sort(names(feature_counts[feature_counts >= min_cohorts]))
      } else if (CONFIG$feature_mode == "intersection") {
        # Pure intersection across cohorts
        common_features_by_tissue[[tissue]] <- sort(Reduce(intersect, tissue_cohorts))
      } else {
        stop("CONFIG$feature_mode must be either 'threshold' or 'intersection'")
      }
    }
    
    log_info(paste("Tissue", tissue, ":", length(common_features_by_tissue[[tissue]]), "common features"))
  }
  
  # Find global intersection across all tissues
  global_features <- Reduce(intersect, common_features_by_tissue)
  log_info(paste("Global intersection:", length(global_features), "features"))
  
  result_obj <- list(
    mode = CONFIG$feature_mode,
    threshold = if (CONFIG$feature_mode == "threshold") CONFIG$feature_threshold else NA,
    by_tissue = common_features_by_tissue,
    global = global_features
  )
  
  save_cache(result_obj, cache_name)  # Cache just the global features
  # saveRDS(result_obj, file.path(CONFIG$output_dir, paste0("common_features_", CONFIG$feature_mode, "_", gsub("\\.", "", as.character(CONFIG$feature_threshold)), ".rds")))
  
  return(global_features)
}

########## Create cell mapping ##########
create_cell_mapping <- function(participants, common_features) {
  cache_name <- "cell_mapping"
  
  if (cache_exists(cache_name) && !CONFIG$force_rerun$cell_mapping) {
    log_info("=== STEP 3: Loading cached cell mapping ===")
    result <- load_cache(cache_name)
    log_info(paste("Loaded", nrow(result), "cells from cache"))
    return(result)
  }
  
  log_info("=== Step 3: Creating cell mapping ===")
  
  cohorts <- names(participants)
  all_cell_data <- list()

  for (cohort in cohorts) {
    log_info(paste("Processing cell mapping for cohort:", cohort))
    
    tryCatch({
      sdrf_file <- file.path(CONFIG$metadata_dir, paste0(cohort, "__", cohort, ".sdrf.txt"))
      cell_meta_file <- file.path(CONFIG$metadata_dir, paste0(cohort, "__", cohort, ".cell_metadata.tsv"))
      barcode_file <- file.path(CONFIG$assay_dir, cohort, paste0(cohort, ".aggregated_filtered_normalised_counts.mtx_cols.gz"))
      
      sdrf <- read_tsv(sdrf_file, col_types = cols(.default = "c"), show_col_types = FALSE)
      cmdf <- read_tsv(cell_meta_file, col_types = cols(.default = "c"), show_col_types = FALSE)
      barcodes <- read_tsv(barcode_file, col_names = "barcode", show_col_types = FALSE)
      
      # Start mapping from barcodes
      mapping <- barcodes %>% mutate(cohort = cohort, original_index = row_number())
      
      # Add cell type from cell_metadata.tsv - our cohorts are OK
      id_col <- grep("^id$|^barcode$", colnames(cmdf), value = TRUE)[1]
      author_cell_type_col <- grep("inferred.*cell.*type.*authors.*labels?", colnames(cmdf), value = TRUE)[1]
      ontology_cell_type_col <- grep("inferred.*cell.*type.*ontology.*labels?", colnames(cmdf), value = TRUE)[1]
      organism_col <- grep("\\borganism\\b", colnames(cmdf), value = TRUE)[1]
      if (cohort == "E-ANND-2") {
        individual_col <- grep("\\bParticipant\\b", colnames(cmdf), value = TRUE)[1]
      } else {
        individual_col <- grep("\\bindividual\\b", colnames(cmdf), value = TRUE)[1] 
      }
      tissue_col <- grep("organism[\\. _-]*part", colnames(cmdf), value = TRUE)[1]
      
      cmdf_sub <- cmdf %>%
        select(barcode = all_of(id_col),
               organism = if(!is.na(organism_col)) all_of(organism_col) else NULL,
               individual = if(!is.na(individual_col)) all_of(individual_col) else NULL,
               tissue_raw = if(!is.na(tissue_col)) all_of(tissue_col) else NULL,
               author_cell_type = if(!is.na(author_cell_type_col)) all_of(author_cell_type_col) else NULL,
               ontology_cell_type = if(!is.na(ontology_cell_type_col)) all_of(ontology_cell_type_col) else NULL) %>%
        distinct()
      
      required_cols <- c("organism", "individual", "tissue_raw", "author_cell_type", "ontology_cell_type")
      for (col in required_cols) {
        if (!col %in% colnames(cmdf_sub)) {
          cmdf_sub[[col]] <- NA_character_
        }
      }
      
      mapping <- mapping %>%
        left_join(cmdf_sub, by = "barcode")
      
      # Filter for human samples and tissues we care about
      cohort_participants <- unlist(participants[[cohort]])
      
      mapping <- mapping %>%
        filter(
          !is.na(organism), organism == "Homo sapiens",
          !is.na(individual), individual %in% cohort_participants,
          !is.na(tissue_raw),
          !is.na(author_cell_type) | !is.na(ontology_cell_type)
        )
      
      mapping$tissue_label <- NA_character_
      
      for (tissue in names(CONFIG$tissues)) {
        if (!tissue %in% names(participants[[cohort]])) next
        
        participant_ids <- participants[[cohort]][[tissue]]
        tissue_values <- tolower(mapping$tissue_raw)
        matches_tissue <- str_detect(tissue_values, tissue)
        matches_participant <- mapping$individual %in% participant_ids
        
        # Apply tissue-specific exclusions if they exist 
        exclusions <- CONFIG$tissue_exclusions[[tissue]]
        if (!is.null(exclusions) && length(exclusions) > 0) {
          exclude_tissue <- str_detect(tissue_values, paste(exclusions, collapse = "|"))
        } else {
          exclude_tissue <- rep(FALSE, nrow(mapping))
        }
        
        include_cells <- matches_tissue & matches_participant & !exclude_tissue
        mapping$tissue_label[include_cells] <- tissue
      }
      
      mapping <- mapping %>%
        mutate(
          cell_type_label = coalesce(ontology_cell_type, author_cell_type) # (prefer ontology labels over author)
        ) %>%
        filter(
          !is.na(tissue_label),
          !is.na(cell_type_label),
          nchar(trimws(cell_type_label)) > 0
        ) %>%
        select(cohort, barcode, original_index, individual, tissue_label, cell_type_label)
      
      # Sanity checks
      if (nrow(mapping) == 0) {
        log_info(paste("No valid cells found for cohort", cohort))
        next
      }
      stopifnot(all(!duplicated(mapping$barcode)))       # each barcode should be unique
      
      # Save cohort-specific mapping
      cohort_file <- file.path(CONFIG$cache_dir, paste0("cell_mapping_", cohort, ".rds"))
      saveRDS(mapping, cohort_file)
      all_cell_data[[cohort]] <- mapping
      
      log_info(paste("Cohort", cohort, ":", nrow(mapping), "cells mapped"))
      
      rm(sdrf, cmdf, barcodes, mapping)
      gc()
      
    }, error = function(e) {
      log_info(paste("Error processing", cohort, ":", e$message), "ERROR")
    })
  }
  
  combined_cell_data <- bind_rows(all_cell_data)
  
  log_info(paste("Total cells BEFORE filtering:", nrow(combined_cell_data)))
  
  # Filter by minimum cells per cell type per tissue
  cell_type_counts <- combined_cell_data %>%
    count(tissue_label, cell_type_label) %>%
    filter(n >= CONFIG$min_cells)
  
  filtered_cell_data <- combined_cell_data %>%
    semi_join(cell_type_counts, by = c("tissue_label", "cell_type_label"))
  
  log_info(paste("Total cells AFTER filtering:", nrow(filtered_cell_data)))
  
  save_cache(filtered_cell_data, cache_name)
  # saveRDS(filtered_cell_data, file.path(CONFIG$output_dir, "cell_mapping.rds"))
  
  return(filtered_cell_data)
}

########## Build expression matrices and create labeled datasets ##########
build_matrices_and_labels <- function(cell_mapping, common_features) {
  log_info("=== Step 4: Building expression matrices and creating labeled datasets ===")
  
  # Create tissue labels
  tissue_labels <- factor(cell_mapping$tissue_label, 
                          levels = names(CONFIG$tissues))
  
  # Create cell type labels  
  cell_type_labels <- factor(cell_mapping$cell_type_label)
  
  # Generate summary statistics
  tissue_summary <- table(tissue_labels)
  cell_type_summary <- table(cell_type_labels)
  log_info(paste("Tissue distribution:", paste(names(tissue_summary), tissue_summary, sep = ":", collapse = ", ")))
  log_info(paste("Total unique cell types:", length(levels(cell_type_labels))))
  
  cohorts <- unique(cell_mapping$cohort)
  all_matrices <- vector("list", length(cohorts))
  names(all_matrices) <- cohorts
  
  for (i in seq_along(cohorts)) {
    cohort <- cohorts[i]
    log_info(paste("Processing cohort", i, "of", length(cohorts), ":", cohort))
    
    tryCatch({
      mtx_file <- file.path(CONFIG$assay_dir, cohort, 
                            paste0(cohort, ".aggregated_filtered_normalised_counts.mtx.gz"))
      row_file <- file.path(CONFIG$assay_dir, cohort, 
                            paste0(cohort, ".aggregated_filtered_normalised_counts.mtx_rows.gz"))
      
      features <- read_tsv(row_file, col_names = "ensembl_id", show_col_types = FALSE)$ensembl_id
      mtx <- readMM(mtx_file) # genes x cells
      log_info(paste("  Raw matrix:", nrow(mtx), "genes x", ncol(mtx), "cells"))
      mtx <- as.csc.matrix(mtx)
      
      # Get cells for this cohort
      cohort_mapping <- cell_mapping %>% filter(cohort == !!cohort)
      barcode_indexes <- cohort_mapping$original_index
      
      # Subset the matrix
      if (CONFIG$feature_mode == "threshold") {
        cohort_mtx <- Matrix(0, nrow = length(common_features), ncol = length(barcode_indexes), sparse = TRUE, 
                             dimnames = list(common_features, cohort_mapping$barcode))
        
        valid_rows_in_cohort_mtx <- which(common_features %in% features)  
        mtx_rows <- match(common_features[valid_rows_in_cohort_mtx], features)  
        mtx_valid <- mtx[mtx_rows, barcode_indexes, drop = FALSE]
        cohort_mtx[valid_rows_in_cohort_mtx, ] <- mtx_valid
        rm(mtx_valid)
      } else if (CONFIG$feature_mode == "intersection") {
        mtx_rows <- match(common_features, features)
        cohort_mtx <- mtx[mtx_rows, barcode_indexes, drop = FALSE]
        rownames(cohort_mtx) <- common_features
        colnames(cohort_mtx) <- cohort_mapping$barcode
      }
      
      rm(mtx)
      gc()
      
      cohort_mtx <- MatrixExtra::t(cohort_mtx) # cells x genes
      log_info(paste("  Final matrix:", nrow(cohort_mtx), "cells x", ncol(cohort_mtx), "features"))

      qs_save(cohort_mtx, file.path(CONFIG$cohort_specific, paste0("MTX_", cohort, "_feature_mode_", CONFIG$feature_mode, 
                                                                   gsub("\\.", "", as.character(CONFIG$feature_threshold)), ".qs")))
      all_matrices[[cohort]] <- cohort_mtx
    
      rm(cohort_mtx, features, cohort_mapping)
      gc()
      
    }, error = function(e) {
      log_info(paste("Error in cohort", cohort, ":", e$message), "ERROR")
    })
  }
  
  combined_matrix <- do.call(MatrixExtra::rbind_csr, all_matrices)
  colnames(combined_matrix) <- common_features
  
  log_info(paste("Combined matrix:", nrow(combined_matrix), "cells x", ncol(combined_matrix), "features"))
  
  log_info("Creating labeled datasets")
  
  feature_names <- colnames(combined_matrix)
  barcode_names <- rownames(combined_matrix)
  
  dataset_name <- generate_dataset_name()
  
  qs_save(tissue_labels, file.path(CONFIG$output_dir, paste0(dataset_name, "_tissue_labels.qs")))
  qs_save(cell_type_labels, file.path(CONFIG$output_dir, paste0(dataset_name, "_celltype_labels.qs")))
  qs_save(feature_names, file.path(CONFIG$output_dir, paste0(dataset_name, "_features.qs")))
  qs_save(barcode_names, file.path(CONFIG$output_dir, paste0(dataset_name, "_barcodes.qs")))
  
  qs_save(combined_matrix, file.path(CONFIG$output_dir, paste0(dataset_name, "_expression_matrix.qs")))
  log_info(paste("Saved expression matrix:", nrow(combined_matrix), "cells x", ncol(combined_matrix), "features"))
  
  dataset_metadata <- list(
    dataset_name = dataset_name,
    n_cells = nrow(combined_matrix),
    n_features = ncol(combined_matrix),
    n_tissues = length(levels(tissue_labels)),
    n_cell_types = length(levels(cell_type_labels)),
    tissue_levels = levels(tissue_labels),
    cell_type_levels = levels(cell_type_labels),
    tissue_distribution = tissue_summary,
    cell_type_distribution = cell_type_summary,
    config = CONFIG,
    created_date = Sys.time(),
    files = list(
      expression_matrix = paste0(dataset_name, "_expression_matrix.qs"),
      tissue_labels = paste0(dataset_name, "_tissue_labels.qs"),
      celltype_labels = paste0(dataset_name, "_celltype_labels.qs"),
      feature_names = paste0(dataset_name, "_features.qs"),
      barcode_names = paste0(dataset_name, "_barcodes.qs")
    )
  )
  
  
  saveRDS(dataset_metadata, file.path(CONFIG$output_dir, paste0(dataset_name, "_metadata.rds")))
  log_info(paste("Final dataset:", dataset_metadata$n_cells, "cells x", dataset_metadata$n_features, "genes"))
  
  rm(combined_matrix, all_matrices)
  gc()
  
  return(dataset_metadata)
}

########## Main pipeline ##########

main <- function() {
  log_info("====================================================================")
  log_info("Starting...")
  log_info(paste("Start time:", Sys.time()))
  log_info(paste("R version:", R.version.string))
  log_info(paste("Feature mode:", CONFIG$feature_mode))
  if (CONFIG$feature_mode == "threshold") {
    log_info(paste("Feature threshold:", CONFIG$feature_threshold))
  }
  log_info(paste("Working directory:", getwd()))
  log_info("====================================================================")
  
  log_info("CACHE STATUS:")
  cache_files <- c("participants", "cell_mapping")
  for (cache_file in cache_files) {
    status <- if (cache_exists(cache_file)) "EXISTS" else "MISSING"
    log_info(paste("-", cache_file, ":", status))
  }
  
  participants <- identify_participants()
  common_features <- find_common_features(participants)
  cell_mapping <- create_cell_mapping(participants, common_features)
  dataset_metadata <- build_matrices_and_labels(cell_mapping, common_features)

  log_info("====================================================================")
  log_info("Completed!")
  log_info(paste("End time:", Sys.time()))
  log_info("====================================================================")
  
  return(dataset_metadata)
}

if (!interactive()) {
  dataset_metadata <- main()
  
  cat("\n====================================================================\n")
  cat("SUMMARY\n")
  cat("====================================================================\n")
  cat("Dataset name:", dataset_metadata$dataset_name, "\n")
  cat("Total cells:", dataset_metadata$n_cells, "\n")
  cat("Total genes:", dataset_metadata$n_features, "\n")
  cat("Feature mode:", CONFIG$feature_mode, "\n")
  cat("Tissues:", dataset_metadata$n_tissues, "\n")
  cat("Cell types:", dataset_metadata$n_cell_types, "\n")
  cat("Output directory:", CONFIG$output_dir, "\n")
  cat("====================================================================\n")
  
  quit(save = "no", status = 0)
}