read_cgm_data <- function(path = here::here("data", "raw", "cgm.csv")) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    dplyr::select(1:3) %>%   # <-- only first 3 columns
    dplyr::rename(
      patient_id = 1,
      glucose = 2,
      time = 3
    ) %>%
    dplyr::mutate(
      patient_id = as.character(patient_id),
      glucose = as.numeric(glucose),
      time = as.numeric(time)
    ) %>%
    dplyr::arrange(patient_id, time)
}

regularize_cgm_grid <- function(data, interval_mins = 5) {
  
  patient_ranges <- data %>%
    dplyr::group_by(patient_id) %>%
    dplyr::summarise(
      min_time = min(time, na.rm = TRUE),
      max_time = max(time, na.rm = TRUE),
      .groups = "drop"
    )
  
  full_grid <- patient_ranges %>%
    dplyr::rowwise() %>%
    dplyr::do(
      tibble::tibble(
        patient_id = .$patient_id,
        time = seq(.$min_time, .$max_time, by = interval_mins)
      )
    ) %>%
    dplyr::ungroup()
  
  full_grid %>%
    dplyr::left_join(data, by = c("patient_id", "time")) %>%
    dplyr::arrange(patient_id, time)
}

label_missing_intervals <- function(data, interval_mins = 5) {
  data %>%
    dplyr::group_by(patient_id) %>%
    dplyr::arrange(time, .by_group = TRUE) %>%
    dplyr::mutate(
      time_diff = dplyr::lead(time) - time,
      gap_steps = ifelse(is.na(time_diff), 0, time_diff / interval_mins),
      missing_start = !is.na(time_diff) & gap_steps > 1,
      miss_duration_mins = ifelse(missing_start, time_diff - interval_mins, NA_real_)
    ) %>%
    dplyr::ungroup()
}

flag_future_missing_onset <- function(data, horizon_steps = 1L) {
  if (horizon_steps != 1L) {
    warning("With gap-based labeling, horizon_steps is being ignored and onset is defined at the current row.")
  }
  
  data %>%
    dplyr::mutate(
      y_onset = as.integer(missing_start)
    )
}