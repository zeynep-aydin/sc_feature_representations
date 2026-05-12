compute_metrics <- function(pred, prob, y_test) {
  conf_matrix <- caret::confusionMatrix(pred, y_test)
  auroc <- tryCatch(
    as.numeric(pROC::multiclass.roc(response = y_test, predictor = prob)$auc),
    error = function(e) NA
  )
  list(
    accuracy = mean(pred == y_test),
    precision = mean(conf_matrix$byClass[, "Precision"], na.rm = TRUE),
    recall = mean(conf_matrix$byClass[, "Recall"], na.rm = TRUE),
    f1_score = mean(conf_matrix$byClass[, "F1"], na.rm = TRUE),
    auroc = auroc,
    confusion_matrix = conf_matrix
  )
}
