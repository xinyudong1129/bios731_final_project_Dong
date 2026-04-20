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
  step_n <- 1L
  
  data |>
    dplyr::group_by(patient_id) |>
    dplyr::mutate(
      is_missing = is.na(glucose),
      miss_run_id = data.table::rleid(is_missing),
      miss_run_length_n = dplyr::if_else(
        is_missing,
        ave(is_missing, miss_run_id, FUN = length),
        0L
      ),
      miss_duration_mins = miss_run_length_n * interval_mins,
      missing_start = is_missing & !dplyr::lag(is_missing, default = FALSE),
      next_missing_start = dplyr::lead(missing_start, default = FALSE)
    ) |>
    dplyr::ungroup()
}

flag_future_missing_onset <- function(data, horizon_steps = 1L) {
  data |>
    dplyr::group_by(patient_id) |>
    dplyr::mutate(
      y_onset = as.integer(dplyr::lead(missing_start, n = horizon_steps, default = FALSE))
    ) |>
    dplyr::ungroup()
}