compute_metrics <- function(pred, prob, y_test) {
  conf_matrix <- caret::confusionMatrix(pred, y_test)
  by_class <- if (is.matrix(conf_matrix$byClass)) conf_matrix$byClass else t(as.matrix(conf_matrix$byClass))
  auroc <- tryCatch(
    as.numeric(pROC::multiclass.roc(response = y_test, predictor = prob)$auc),
    error = function(e) NA
  )
  list(
    accuracy = mean(pred == y_test),
    precision = mean(by_class[, "Precision"], na.rm = TRUE),
    recall = mean(by_class[, "Recall"], na.rm = TRUE),
    f1_score = mean(by_class[, "F1"], na.rm = TRUE),
    auroc = auroc,
    confusion_matrix = conf_matrix
  )
}
