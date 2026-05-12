train_glmnet <- function(X_train, y_train) {
  t_start <- Sys.time()
  model <- glmnet::glmnet(x = X_train, y = y_train, family = "multinomial", alpha = 1)
  model_time <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  best_lambda <- model$lambda[which.max(model$dev.ratio)]
  list(
    model = model,
    best_lambda = best_lambda,
    model_time = model_time,
    model_params = list(
      family = "multinomial",
      alpha = 1,
      best_lambda = best_lambda,
      lambda_selection = "max_dev_ratio"
    )
  )
}

train_svm <- function(X_train, y_train) {
  t_start <- Sys.time()
  model <- e1071::svm(x = X_train, y = y_train, kernel = "linear",
                      type = "C-classification", scale = TRUE, probability = TRUE)
  model_time <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  list(
    model = model,
    model_time = model_time,
    model_params = list(
      kernel = "linear",
      type = "C-classification",
      scale = TRUE,
      probability = TRUE
    )
  )
}

predict_glmnet <- function(model, best_lambda, X_test, y_test) {
  t_start <- Sys.time()
  pred <- factor(
    predict(model, newx = X_test, s = best_lambda, type = "class"),
    levels = levels(y_test)
  )
  prob <- predict(model, newx = X_test, s = best_lambda, type = "response")[, , 1]
  list(pred = pred, prob = prob, predict_time = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
}

predict_svm <- function(model, X_test, y_test) {
  t_start <- Sys.time()
  pred <- predict(model, newdata = as.matrix(X_test))
  prob_raw <- predict(model, newdata = as.matrix(X_test), probability = TRUE)
  list(pred = pred, prob = attr(prob_raw, "probabilities"), predict_time = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
}
