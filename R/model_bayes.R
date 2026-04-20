fit_bayes_mcmc_onset <- function(train_data,
                                 chains = 4,
                                 iter = 2000,
                                 seed = 123) {
  stopifnot(requireNamespace("rstanarm", quietly = TRUE))
  
  rstanarm::stan_glmer(
    y_onset ~ mean_glucose + sd_glucose + slope_glucose +
      delta_last_first + pct_high + pct_low +
      hour_of_day + day_of_week + (1 | patient_id),
    data = train_data,
    family = binomial(link = "logit"),
    chains = chains,
    iter = iter,
    seed = seed,
    refresh = 0
  )
}

fit_bayes_mcmc_duration <- function(train_data,
                                    chains = 4,
                                    iter = 2000,
                                    seed = 123) {
  stopifnot(requireNamespace("rstanarm", quietly = TRUE))
  
  train_pos <- train_data |>
    dplyr::filter(y_onset == 1, !is.na(y_duration)) |>
    dplyr::mutate(log_duration = log1p(y_duration))
  
  rstanarm::stan_lmer(
    log_duration ~ mean_glucose + sd_glucose + slope_glucose +
      delta_last_first + pct_high + pct_low +
      hour_of_day + day_of_week + (1 | patient_id),
    data = train_pos,
    chains = chains,
    iter = iter,
    seed = seed,
    refresh = 0
  )
}

fit_bayes_vb_onset <- function(train_data, seed = 123) {
  stopifnot(requireNamespace("rstanarm", quietly = TRUE))
  
  rstanarm::stan_glmer(
    y_onset ~ mean_glucose + sd_glucose + slope_glucose +
      delta_last_first + pct_high + pct_low +
      hour_of_day + day_of_week + (1 | patient_id),
    data = train_data,
    family = binomial(link = "logit"),
    algorithm = "meanfield",
    seed = seed,
    refresh = 0
  )
}

fit_bayes_vb_duration <- function(train_data, seed = 123) {
  stopifnot(requireNamespace("rstanarm", quietly = TRUE))
  
  train_pos <- train_data |>
    dplyr::filter(y_onset == 1, !is.na(y_duration)) |>
    dplyr::mutate(log_duration = log1p(y_duration))
  
  rstanarm::stan_lmer(
    log_duration ~ mean_glucose + sd_glucose + slope_glucose +
      delta_last_first + pct_high + pct_low +
      hour_of_day + day_of_week + (1 | patient_id),
    data = train_pos,
    algorithm = "meanfield",
    seed = seed,
    refresh = 0
  )
}

predict_bayes_onset <- function(fit, new_data) {
  as.numeric(stats::predict(fit, newdata = new_data, type = "response"))
}

predict_bayes_duration <- function(fit, new_data) {
  pred <- as.numeric(stats::predict(fit, newdata = new_data))
  pmax(0, exp(pred) - 1)
}

posterior_interval_duration <- function(fit, new_data, prob = 0.95) {
  draws <- posterior_epred_safe(fit, new_data)
  alpha <- (1 - prob) / 2
  qs <- apply(draws, 2, stats::quantile, probs = c(alpha, 1 - alpha))
  tibble::tibble(
    lower = pmax(0, exp(qs[1, ]) - 1),
    upper = pmax(0, exp(qs[2, ]) - 1)
  )
}

posterior_epred_safe <- function(fit, new_data) {
  if ("stanreg" %in% class(fit)) {
    rstanarm::posterior_epred(fit, newdata = new_data)
  } else {
    stop("Unsupported fit object.")
  }
}