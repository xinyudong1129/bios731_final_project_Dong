fit_glm_em_onset <- function(train_data, k_components = 2) {
  if (length(unique(stats::na.omit(train_data$y_onset))) < 2) {
    stop("y_onset must contain both 0 and 1 in the training data.")
  }
  
  stats::glm(
    y_onset ~ mean_glucose + sd_glucose + slope_glucose +
      delta_last_first + pct_high + pct_low,
    data = train_data,
    family = stats::binomial()
  )
}

predict_glm_em_onset <- function(fit, new_data) {
  as.numeric(stats::predict(fit, newdata = new_data, type = "response"))
}

fit_glm_em_duration <- function(train_data) {
  train_pos <- train_data %>%
    dplyr::filter(y_onset == 1, !is.na(y_duration)) %>%
    dplyr::mutate(log_duration = log1p(y_duration))
  
  if (nrow(train_pos) < 10) {
    return(NULL)
  }
  
  stats::lm(
    log_duration ~ mean_glucose + sd_glucose + slope_glucose +
      delta_last_first + pct_high + pct_low,
    data = train_pos
  )
}

predict_glm_em_duration <- function(fit, new_data) {
  if (is.null(fit)) {
    return(rep(NA_real_, nrow(new_data)))
  }
  
  pred <- stats::predict(fit, newdata = new_data)
  pmax(0, exp(pred) - 1)
}