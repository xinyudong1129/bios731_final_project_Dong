fit_bayes_mcmc_onset <- function(train_data,
                                 chains = 2,
                                 iter = 1000,
                                 seed = 123) {
  stopifnot(requireNamespace("rstanarm", quietly = TRUE))
  
  rstanarm::stan_glm(
    y_onset ~ mean_glucose + sd_glucose + slope_glucose +
      delta_last_first + pct_high + pct_low,
    data = train_data,
    family = binomial(link = "logit"),
    chains = chains,
    iter = iter,
    seed = seed,
    refresh = 0
  )
}

predict_bayes_onset <- function(fit, new_data) {
  draws <- rstanarm::posterior_epred(fit, newdata = new_data)
  colMeans(draws)
}

fit_bayes_mcmc_duration <- function(train_data,
                                    chains = 2,
                                    iter = 1000,
                                    seed = 123) {
  stopifnot(requireNamespace("rstanarm", quietly = TRUE))
  
  train_pos <- train_data %>%
    dplyr::filter(y_onset == 1, !is.na(y_duration)) %>%
    dplyr::mutate(log_duration = log1p(y_duration))
  
  if (nrow(train_pos) < 30) {
    message("Skipping Bayesian duration: too few rows")
    return(NULL)
  }
  
  rstanarm::stan_glm(
    log_duration ~ mean_glucose + sd_glucose + slope_glucose +
      delta_last_first + pct_high + pct_low,
    data = train_pos,
    family = gaussian(),
    chains = chains,
    iter = iter,
    seed = seed,
    refresh = 0
  )
}

predict_bayes_duration <- function(fit, new_data) {
  if (is.null(fit)) {
    return(rep(NA_real_, nrow(new_data)))
  }
  
  draws <- rstanarm::posterior_epred(fit, newdata = new_data)
  pred_log <- colMeans(draws)
  pmax(0, exp(pred_log) - 1)
}

posterior_interval_duration <- function(fit, new_data, prob = 0.95) {
  if (is.null(fit)) {
    return(tibble::tibble(
      lower = rep(NA_real_, nrow(new_data)),
      upper = rep(NA_real_, nrow(new_data))
    ))
  }
  
  draws <- rstanarm::posterior_epred(fit, newdata = new_data)
  alpha <- (1 - prob) / 2
  
  qs <- apply(draws, 2, stats::quantile, probs = c(alpha, 1 - alpha), na.rm = TRUE)
  
  tibble::tibble(
    lower = pmax(0, exp(qs[1, ]) - 1),
    upper = pmax(0, exp(qs[2, ]) - 1)
  )
}