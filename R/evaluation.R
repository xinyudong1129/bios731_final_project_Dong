compute_auroc <- function(truth, score) {
  keep <- !is.na(truth) & !is.na(score)
  truth <- truth[keep]
  score <- score[keep]
  
  # AUROC is undefined if only one class is present
  if (length(unique(truth)) < 2) {
    return(NA_real_)
  }
  
  as.numeric(
    pROC::roc(
      response = truth,
      predictor = score,
      quiet = TRUE
    )$auc
  )
}

compute_auprc <- function(truth, score) {
  keep <- !is.na(truth) & !is.na(score)
  truth <- truth[keep]
  score <- score[keep]
  
  # AUPRC is undefined if one class is missing
  if (length(unique(truth)) < 2) {
    return(NA_real_)
  }
  
  pos <- score[truth == 1]
  neg <- score[truth == 0]
  
  if (length(pos) == 0 || length(neg) == 0) {
    return(NA_real_)
  }
  
  PRROC::pr.curve(
    scores.class0 = pos,
    scores.class1 = neg
  )$auc.integral
}

compute_rmse <- function(truth, pred) {
  keep <- !is.na(truth) & !is.na(pred)
  
  if (sum(keep) == 0) {
    return(NA_real_)
  }
  
  sqrt(mean((truth[keep] - pred[keep])^2))
}

compute_interval_coverage <- function(truth, lower, upper) {
  keep <- !is.na(truth) & !is.na(lower) & !is.na(upper)
  
  if (sum(keep) == 0) {
    return(NA_real_)
  }
  
  mean(truth[keep] >= lower[keep] & truth[keep] <= upper[keep])
}

evaluate_predictions <- function(data,
                                 onset_prob,
                                 duration_pred,
                                 duration_lower = NULL,
                                 duration_upper = NULL,
                                 method_name = "unknown") {
  if (length(unique(data$y_onset)) < 2) {
    return(tibble::tibble(
      method = method_name,
      auroc = NA_real_,
      auprc = NA_real_,
      rmse_duration = NA_real_,
      coverage_95 = NA_real_
    ))
  }
  # Onset metrics
  auroc_val <- compute_auroc(
    truth = data$y_onset,
    score = onset_prob
  )
  
  auprc_val <- compute_auprc(
    truth = data$y_onset,
    score = onset_prob
  )
  
  # Duration metrics only among true onset rows with observed duration
  dur_idx <- data$y_onset == 1 & !is.na(data$y_duration)
  
  rmse_val <- compute_rmse(
    truth = data$y_duration[dur_idx],
    pred = duration_pred[dur_idx]
  )
  
  out <- tibble::tibble(
    method = method_name,
    auroc = auroc_val,
    auprc = auprc_val,
    rmse_duration = rmse_val,
    coverage_95 = NA_real_
  )
  
  # Coverage only if interval bounds are supplied
  if (!is.null(duration_lower) && !is.null(duration_upper)) {
    out$coverage_95 <- compute_interval_coverage(
      truth = data$y_duration[dur_idx],
      lower = duration_lower[dur_idx],
      upper = duration_upper[dur_idx]
    )
  }
  
  out
}