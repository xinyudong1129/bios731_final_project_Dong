calc_window_features <- function(x, timestamps) {
  n <- length(x)
  idx <- seq_len(n)
  
  if (all(is.na(x))) {
    return(tibble::tibble(
      mean_glucose = NA_real_,
      sd_glucose = NA_real_,
      min_glucose = NA_real_,
      max_glucose = NA_real_,
      range_glucose = NA_real_,
      slope_glucose = NA_real_,
      last_glucose = NA_real_,
      delta_last_first = NA_real_,
      pct_high = NA_real_,
      pct_low = NA_real_,
      hour_of_day = NA_real_,
      day_of_week = NA_real_
    ))
  }
  
  fit <- stats::lm(x ~ idx)
  last_time <- max(timestamps, na.rm = TRUE)
  
  tibble::tibble(
    mean_glucose = mean(x, na.rm = TRUE),
    sd_glucose = stats::sd(x, na.rm = TRUE),
    min_glucose = min(x, na.rm = TRUE),
    max_glucose = max(x, na.rm = TRUE),
    range_glucose = max(x, na.rm = TRUE) - min(x, na.rm = TRUE),
    slope_glucose = unname(stats::coef(fit)[2]),
    last_glucose = x[length(x)],
    delta_last_first = x[length(x)] - x[1],
    pct_high = mean(x > 180, na.rm = TRUE),
    pct_low = mean(x < 70, na.rm = TRUE),
    hour_of_day = floor((last_time %% (24 * 60)) / 60),
    day_of_week = floor(last_time / (24 * 60)) %% 7 + 1
  )
}
make_pre_missing_windows <- function(data,
                                     window_size = 12L,
                                     horizon_steps = 1L,
                                     interval_mins = 5) {
  data <- data |>
    flag_future_missing_onset(horizon_steps = horizon_steps)
  
  out <- data |>
    dplyr::group_by(patient_id) |>
    dplyr::group_modify(function(df, key) {
      n <- nrow(df)
      rows <- vector("list", n)
      
      for (i in seq(window_size, n - horizon_steps)) {
        window_idx <- (i - window_size + 1):i
        future_idx <- i + horizon_steps
        
        x <- df$glucose[window_idx]
        tt <- df$time[window_idx]
        
        if (any(is.na(x))) next
        
        feats <- calc_window_features(x, tt)
        
        rows[[i]] <- dplyr::bind_cols(
          tibble::tibble(
            end_time = df$time[i],
            y_onset = df$y_onset[i],
            y_duration = ifelse(df$y_onset[i] == 1,
                                df$miss_duration_mins[future_idx],
                                NA_real_)
          ),
          feats
        )
      }
      
      dplyr::bind_rows(rows)
    }) |>
    dplyr::ungroup()
  
  out
}
build_analysis_dataset <- function(raw_data,
                                   interval_mins = 5,
                                   window_size = 12L,
                                   horizon_steps = 1L) {
  raw_data |>
    regularize_cgm_grid(interval_mins = interval_mins) |>
    label_missing_intervals(interval_mins = interval_mins) |>
    make_pre_missing_windows(
      window_size = window_size,
      horizon_steps = horizon_steps,
      interval_mins = interval_mins
    ) |>
    dplyr::mutate(
      hour_of_day = factor(hour_of_day),
      day_of_week = factor(day_of_week)
    )
}