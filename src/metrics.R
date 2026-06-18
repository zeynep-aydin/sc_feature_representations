# Macro precision/recall/F1 over test-present classes (zero_division = 0)
# accuracy is micro (NA preds counted as wrong)
# AUROC over test-present classes only
compute_metrics <- function(pred, prob, y_test) {
  # Align pred and truth to a common level set so the confusion matrix is square and name-aligned
  # pred can contain trained classes that are absent from the test set
  lev <- union(levels(factor(y_test)), levels(factor(pred)))
  predf <- factor(pred, levels = lev)
  truth <- factor(y_test, levels = lev)
  cm <- table(pred = predf, truth = truth)
  tp <- diag(cm)
  support <- colSums(cm)                 # test instances per class
  precision <- tp / rowSums(cm)          # NaN when class never predicted
  recall <- tp / support                 # NaN when class absent from test
  precision[is.nan(precision)] <- 0      # zero_division = 0
  recall[is.nan(recall)] <- 0
  f1 <- 2 * precision * recall / (precision + recall)
  f1[is.nan(f1)] <- 0
  present <- support > 0                 # macro-average over test-present classes
  # AUROC over test-present classes: match prob columns to the response levels.
  # make_split guarantees every test class is trained, so prob always has these columns.
  auroc <- tryCatch({
    truth_present <- droplevels(truth)
    prob_scored <- prob[, levels(truth_present), drop = FALSE]
    as.numeric(pROC::multiclass.roc(response = truth_present, predictor = prob_scored)$auc)
  }, error = function(e) NA)
  list(
    accuracy = mean(!is.na(pred) & as.character(pred) == as.character(y_test)),
    precision = mean(precision[present]),
    recall = mean(recall[present]),
    f1_score = mean(f1[present]),
    auroc = auroc,
    confusion_matrix = cm
  )
}
