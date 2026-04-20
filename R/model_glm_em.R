#This baseline uses an EM-based latent-class logistic model via flexmix for onset, 
# plus a simple interpretable regression for log-duration.


fit_glm_em_onset <- function(train_data, k_components = 2) {
  stopifnot(requireNamespace("flexmix", quietly = TRUE))
  
  form <- stats::as.formula(
    y_onset ~ mean_glucose + sd_glucose + slope_glucose +
      delta_last_first + pct_high + pct_low +
      hour_of_day + day_of_week
  )
  
  flexmix::flexmix(
    formula = form,
    data = train_data,
    k = k_components,
    model = flexmix::FLXMRglm(family = "binomial")
  )
}

predict_glm_em_onset <- function(fit, new_data) {
  post <- flexmix::posterior(fit, newdata = new_data)
  comp_probs <- sapply(seq_len(ncol(post)), function(j) {
    stats::predict(fit@components[[j]], newdata = new_data, type = "response")
  })
  rowSums(post * comp_probs)
}

fit_glm_em_duration <- function(train_data) {
  train_pos <- train_data |>
    dplyr::filter(y_onset == 1, !is.na(y_duration)) |>
    dplyr::mutate(log_duration = log1p(y_duration))
  
  stats::lm(
    log_duration ~ mean_glucose + sd_glucose + slope_glucose +
      delta_last_first + pct_high + pct_low +
      hour_of_day + day_of_week,
    data = train_pos
  )
}

predict_glm_em_duration <- function(fit, new_data) {
  pred <- stats::predict(fit, newdata = new_data)
  pmax(0, exp(pred) - 1)
}