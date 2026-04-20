read_cgm_data <- function(path = here::here("data", "raw", "cgm.csv")) {
  readr::read_csv(
    path,
    col_types = readr::cols(
      patient_id = readr::col_character(),
      time = readr::col_datetime(),
      glucose = readr::col_double()
    )
  ) |>
    dplyr::arrange(patient_id, time)
}

regularize_cgm_grid <- function(data, interval_mins = 5) {
  interval_sec <- interval_mins * 60
  
  data |>
    dplyr::group_by(patient_id) |>
    dplyr::summarise(
      time = seq(min(time, na.rm = TRUE),
                 max(time, na.rm = TRUE),
                 by = interval_sec),
      .groups = "drop"
    ) |>
    dplyr::left_join(data, by = c("patient_id", "time")) |>
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