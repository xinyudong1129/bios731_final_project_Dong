compute_auroc <- function(truth, score) {
  pROC::roc(response = truth, predictor = score, quiet = TRUE)$auc |>
    as.numeric()
}

compute_auprc <- function(truth, score) {
  pos <- score[truth == 1]
  neg <- score[truth == 0]
  PRROC::pr.curve(scores.class0 = pos, scores.class1 = neg)$auc.integral
}

compute_rmse <- function(truth, pred) {
  sqrt(mean((truth - pred)^2, na.rm = TRUE))
}

compute_interval_coverage <- function(truth, lower, upper) {
  mean(truth >= lower & truth <= upper, na.rm = TRUE)
}

evaluate_predictions <- function(data,
                                 onset_prob,
                                 duration_pred,
                                 duration_lower = NULL,
                                 duration_upper = NULL,
                                 method_name = "unknown") {
  out <- tibble::tibble(
    method = method_name,
    auroc = compute_auroc(data$y_onset, onset_prob),
    auprc = compute_auprc(data$y_onset, onset_prob),
    rmse_duration = compute_rmse(
      truth = data$y_duration[data$y_onset == 1],
      pred = duration_pred[data$y_onset == 1]
    )
  )
  
  if (!is.null(duration_lower) && !is.null(duration_upper)) {
    out$coverage_95 <- compute_interval_coverage(
      truth = data$y_duration[data$y_onset == 1],
      lower = duration_lower[data$y_onset == 1],
      upper = duration_upper[data$y_onset == 1]
    )
  }
  
  out
}