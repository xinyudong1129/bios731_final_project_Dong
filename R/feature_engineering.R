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
      pct_low = NA_real_
    ))
  }
  
  # simple linear interpolation inside the window if there are a few NAs
  if (any(is.na(x))) {
    nonmiss <- which(!is.na(x))
    
    # cannot interpolate if too few observed points
    if (length(nonmiss) < 2) {
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
        pct_low = NA_real_
      ))
    }
    
    x <- stats::approx(
      x = nonmiss,
      y = x[nonmiss],
      xout = seq_along(x),
      method = "linear",
      rule = 2
    )$y
  }
  
  fit <- stats::lm(x ~ idx)
  
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
    pct_low = mean(x < 70, na.rm = TRUE)
  )
}


make_pre_missing_windows <- function(data,
                                     window_size = 6L,
                                     horizon_steps = 1L,
                                     interval_mins = 5,
                                     max_missing_in_window = 2L,
                                     time_tolerance = 1) {
  data <- data %>%
    flag_future_missing_onset(horizon_steps = horizon_steps)
  
  empty_out <- tibble::tibble(
    end_time = numeric(),
    y_onset = integer(),
    y_duration = numeric(),
    mean_glucose = numeric(),
    sd_glucose = numeric(),
    min_glucose = numeric(),
    max_glucose = numeric(),
    range_glucose = numeric(),
    slope_glucose = numeric(),
    last_glucose = numeric(),
    delta_last_first = numeric(),
    pct_high = numeric(),
    pct_low = numeric()
  )
  
  out <- data %>%
    dplyr::group_by(patient_id) %>%
    dplyr::group_modify(function(df, key) {
      n <- nrow(df)
      rows <- vector("list", n)
      
      for (i in seq(window_size, n)) {
        window_idx <- (i - window_size + 1):i
        x <- df$glucose[window_idx]
        tt <- df$time[window_idx]
        
        # allow a small amount of missingness in the window
        if (sum(is.na(x)) > max_missing_in_window) next
        
        # allow slight irregularity in spacing
        #if (length(tt) > 1 && any(abs(diff(tt) - interval_mins) > time_tolerance)) next
        
        feats <- calc_window_features(x, tt)
        
        # skip if feature extraction still failed
        if (all(is.na(feats))) next
        
        rows[[i]] <- dplyr::bind_cols(
          tibble::tibble(
            end_time = df$time[i],
            y_onset = df$y_onset[i],
            y_duration = ifelse(df$y_onset[i] == 1,
                                df$miss_duration_mins[i],
                                NA_real_)
          ),
          feats
        )
      }
      
      rows <- rows[!vapply(rows, is.null, logical(1))]
      if (length(rows) == 0) return(empty_out)
      
      dplyr::bind_rows(rows)
    }) %>%
    dplyr::ungroup()
  
  out
}


build_analysis_dataset <- function(raw_data,
                                   interval_mins = 5,
                                   window_size = 6L,
                                   horizon_steps = 1L,
                                   max_missing_in_window = 2L,
                                   time_tolerance = 1) {
  out <- raw_data %>%
    label_missing_intervals(interval_mins = interval_mins) %>%
    make_pre_missing_windows(
      window_size = window_size,
      horizon_steps = horizon_steps,
      interval_mins = interval_mins,
      max_missing_in_window = max_missing_in_window,
      time_tolerance = time_tolerance
    )
  
  if (nrow(out) == 0) {
    stop("No valid complete windows were created.")
  }
  
  out
}