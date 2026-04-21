simulate_cgm_panel <- function(n_patients = 50,
                               n_time = 288,
                               interval_mins = 5,
                               seed = 123) {
  set.seed(seed)
  
  patient_ids <- sprintf("P%03d", seq_len(n_patients))
  start_time <- 0
  
  sim_list <- lapply(patient_ids, function(id) {
    t_idx <- seq_len(n_time)
    tt <- start_time + (t_idx - 1) * interval_mins
    
    patient_intercept <- rnorm(1, 110, 15)
    circadian <- 15 * sin(2 * pi * t_idx / (24 * 60 / interval_mins))
    ar_noise <- as.numeric(stats::arima.sim(model = list(ar = 0.7), n = n_time, sd = 8))
    glucose_true <- patient_intercept + circadian + ar_noise
    
    lag_glucose <- dplyr::lag(glucose_true, default = glucose_true[1])
    
    # less aggressive onset probability
    p_missing_start <- plogis(-6 + 0.01 * lag_glucose)
    missing_start <- rbinom(n_time, 1, p_missing_start)
    
    # shorter durations
    duration_steps <- ifelse(
      missing_start == 1,
      sample(2:6, n_time, replace = TRUE),
      0
    )
    
    glucose_obs <- glucose_true
    i <- 1
    while (i <= n_time) {
      if (missing_start[i] == 1) {
        end_i <- min(n_time, i + duration_steps[i] - 1)
        glucose_obs[i:end_i] <- NA_real_
        i <- end_i + 1
      } else {
        i <- i + 1
      }
    }
    
    tibble::tibble(
      patient_id = id,
      time = tt,
      glucose = glucose_obs
    ) %>%
      dplyr::filter(!is.na(glucose))
  })
  
  dplyr::bind_rows(sim_list)
}