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

# donor-holdout split with nested cell-level training fractions.
# reserve a donor-pure test set first, independent of train_frac
# the remaining donors form the training pool
make_split <- function(data, labels, batch_labels, train_frac, run_id,
                       train_fracs = c(0.40, 0.60, 0.80)) {
  if (length(unique(batch_labels)) < 2)
    stop("make_split expects >= 2 donors; got a single-donor dataset")
  set.seed(2969 + run_id)
  max_frac <- max(train_fracs)
  labels_chr <- as.character(labels)
  total <- length(labels)

  # prefix length whose cumulative fraction is closest to target (cum_frac monotone increasing)
  prefix_n <- function(cum_frac, target) which.min(abs(cum_frac - target))

  # donor-pure test set = donors beyond the max-fraction prefix.
  don_sizes <- table(batch_labels)
  ordered <- sample(names(don_sizes))                        # seeded random donor order
  cum_frac <- cumsum(as.integer(don_sizes[ordered])) / total
  n_pool <- prefix_n(cum_frac, max_frac)
  pool_don <- ordered[seq_len(n_pool)]
  test_don <- setdiff(ordered, pool_don)
  test_idx <- which(batch_labels %in% test_don)

  # nested cell-level training prefixes drawn from the pool (capped by pool size).
  pool_idx <- sample(which(batch_labels %in% pool_don))       # seeded cell shuffle
  n_train <- min(round(train_frac * total), length(pool_idx))
  train_idx <- pool_idx[seq_len(n_train)]

  log_info(sprintf("Donor-holdout split (run %d, frac %.2f): train=%d cells (%.1f%%), test=[%s] (%.1f%%)",
                   run_id, train_frac, length(train_idx), 100 * length(train_idx) / total,
                   paste(sort(test_don), collapse = "+"), 100 * length(test_idx) / total))

  # evaluable classes: >= 2 training cells in the smallest-fraction prefix
  # identical across train_fracs for a run_id
  min_frac <- min(train_fracs)
  n_eval <- min(round(min_frac * total), length(pool_idx))
  evaltrain_idx <- pool_idx[seq_len(n_eval)]
  eval_counts <- table(labels_chr[evaltrain_idx])
  eval_classes <- names(eval_counts)[eval_counts >= 2]

  # restrict the returned indices to the evaluable classes so test-only/under-supported classes never reach the model
  actual_train_fraction <- length(train_idx) / total
  dropped <- setdiff(unique(labels_chr[c(train_idx, test_idx)]), eval_classes)
  if (length(dropped) > 0)
    log_info(sprintf("Dropping %d non-eval class(es): %s",
                     length(dropped), paste(sort(dropped), collapse = ", ")))
  train_idx <- train_idx[labels_chr[train_idx] %in% eval_classes]
  test_idx <- test_idx[labels_chr[test_idx]  %in% eval_classes]

  # every test-present class must be realised in the returned train set
  stopifnot(all(labels_chr[test_idx] %in% labels_chr[train_idx]))

  list(
    input_cells = nrow(data),
    input_genes = ncol(data),
    train_indices = train_idx,
    test_indices = test_idx,
    eval_classes = eval_classes,
    donor_order = ordered,
    test_donors = test_don,
    actual_train_fraction = actual_train_fraction,
    eval_train_fraction = length(evaltrain_idx) / total,
    split_seed = 2969 + run_id,
    run_id = run_id
  )
}
