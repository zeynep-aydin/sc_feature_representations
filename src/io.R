load_dataset <- function(dataset, data_dir, feature_mode, task) {
  if (dataset == "scea") {
    mtx <- qs2::qs_read(file.path(data_dir, paste0("scea_8tissues_", feature_mode, "_expression_matrix.qs")))
    labels <- qs2::qs_read(file.path(data_dir, paste0("scea_8tissues_", feature_mode, "_", task, "_labels.qs")))
    batch_labels <- qs2::qs_read(file.path(data_dir, paste0("scea_8tissues_", feature_mode, "_cohort_labels.qs")))
  } else {
    mtx <- qs2::qs_read(file.path(data_dir, paste0(dataset, "_intersection_expression_matrix.qs")))
    labels <- qs2::qs_read(file.path(data_dir, paste0(dataset, "_intersection_", task, "_labels.qs")))
    batch_labels <- qs2::qs_read(file.path(data_dir, paste0(dataset, "_intersection_batch_labels.qs")))
  }
  list(mtx = mtx, labels = labels, batch_labels = batch_labels)
}

make_split <- function(data, labels, batch_labels, dataset, train_frac, run_id) {
  set.seed(2969 + run_id)
  unique_groups <- unique(batch_labels)

  if (dataset %in% names(FIXED_SPLITS)) {
    split_key <- paste0(round(train_frac * 100), "%")
    entry_name <- paste0("S", sprintf("%02d", run_id))
    entry <- FIXED_SPLITS[[dataset]][[split_key]][[entry_name]]
    if (is.null(entry)) {
      stop(sprintf("No fixed split found for dataset=%s, train_pct=%s, run_id=%d",
                   dataset, split_key, run_id))
    }
    train_idx <- which(batch_labels %in% entry$train)
    test_idx <- which(batch_labels %in% entry$test)
    log_info(sprintf("Using fixed split %s (%s): %d train / %d test cells",
                     entry_name, split_key, length(train_idx), length(test_idx)))

  } else if (length(unique_groups) == 1) {
    shuffle_idx <- sample(seq_along(labels))
    train_idx <- shuffle_idx[1:floor(train_frac * length(labels))]
    test_idx <- shuffle_idx[(length(train_idx) + 1):length(labels)]

  } else {
    stop(sprintf(
      "Dataset '%s' has multiple groups but no fixed splits defined. Add to fixed_splits.R.",
      dataset
    ))
  }

  list(
    input_cells = nrow(data),
    input_genes = ncol(data),
    train_indices = train_idx,
    test_indices = test_idx,
    split_seed = 2969 + run_id,
    run_id = run_id
  )
}
