make_project_dirs <- function() {
  dirs <- c(
    here::here("data", "raw"),
    here::here("data", "processed"),
    here::here("results", "metrics"),
    here::here("results", "predictions"),
    here::here("results", "models"),
    here::here("results", "simulation"),
    here::here("outputs", "figures"),
    here::here("outputs", "tables")
  )
  
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

save_rds_here <- function(object, ...) {
  path <- here::here(...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, path)
  invisible(path)
}

read_rds_here <- function(...) {
  readRDS(here::here(...))
}

write_csv_here <- function(data, ...) {
  path <- here::here(...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(data, path)
  invisible(path)
}