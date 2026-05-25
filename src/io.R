load_dataset <- function(dataset, data_dir, task = "celltype", label_level = NULL) {
  mtx <- qs2::qs_read(file.path(data_dir, "counts.qs"))
  label_file <- if (task == "tissue") {
    "tissue.qs"
  } else if (!is.null(label_level)) {
    paste0("labels_", label_level, ".qs")
  } else {
    "labels.qs"
  }
  labels <- qs2::qs_read(file.path(data_dir, label_file))
  donors <- qs2::qs_read(file.path(data_dir, "donors.qs"))
  list(mtx = mtx, labels = labels, batch_labels = donors)
}

make_split <- function(data, labels, batch_labels, train_frac, run_id) {
  set.seed(2969 + run_id)
  unique_groups <- unique(batch_labels)

  if (length(unique_groups) == 1) {
    shuffle_idx <- sample(seq_along(labels))
    train_idx <- shuffle_idx[seq_len(floor(train_frac * length(labels)))]
    test_idx <- shuffle_idx[(length(train_idx) + 1):length(labels)]
  } else {
    don_sizes <- table(batch_labels)
    total <- length(batch_labels)
    n <- length(unique_groups)

    if (n <= 20) {
      # Exhaustive: within frac_tol of target, minimize test-only classes then fraction error.
      frac_tol <- 0.01
      don_type_map <- tapply(as.character(labels), batch_labels, unique)
      best_cov_score <- c(Inf, Inf)
      best_train_don <- NULL
      for (mask in seq_len(2^n - 2)) {
        bits <- as.logical(intToBits(mask)[seq_len(n)])
        train_don <- unique_groups[bits]
        frac_err <- abs(sum(don_sizes[train_don]) / total - train_frac)
        if (frac_err <= frac_tol) {
          test_don <- unique_groups[!bits]
          train_types <- unique(unlist(don_type_map[train_don]))
          test_types <- unique(unlist(don_type_map[test_don]))
          n_test_only <- length(setdiff(test_types, train_types))
          score <- c(n_test_only, frac_err)
          if (score[1] < best_cov_score[1] ||
              (score[1] == best_cov_score[1] && score[2] < best_cov_score[2])) {
            best_cov_score <- score
            best_train_don <- train_don
          }
        }
      }
      if (is.null(best_train_don))
        stop(sprintf("No donor partition within %.0fpp of target fraction %.0f%%",
                     frac_tol * 100, train_frac * 100))
    } else {
      # Greedy: seed-based random donor permutation, assign smallest donors to test
      don_df <- data.frame(don = names(don_sizes), n = as.integer(don_sizes),
                           stringsAsFactors = FALSE)
      don_df <- don_df[sample(nrow(don_df)), ]
      don_df <- don_df[order(don_df$n), ]
      target_test <- (1 - train_frac) * total
      test_don <- character(0)
      test_cells <- 0
      for (i in seq_len(nrow(don_df))) {
        if (test_cells >= target_test) break
        test_don <- c(test_don, don_df$don[i])
        test_cells <- test_cells + don_df$n[i]
      }
      best_train_don <- setdiff(names(don_sizes), test_don)
    }

    train_idx <- which(batch_labels %in% best_train_don)
    test_idx <- which(!batch_labels %in% best_train_don)
    log_info(sprintf("Donor split: train=[%s] test=[%s] achieved=%.1f%%",
                     paste(sort(best_train_don), collapse = "+"),
                     paste(sort(setdiff(unique_groups, best_train_don)), collapse = "+"),
                     100 * length(train_idx) / total))
  }

  list(
    input_cells = nrow(data),
    input_genes = ncol(data),
    train_indices = train_idx,
    test_indices = test_idx,
    actual_train_fraction = length(train_idx) / length(labels),
    split_seed = 2969 + run_id,
    run_id = run_id
  )
}
