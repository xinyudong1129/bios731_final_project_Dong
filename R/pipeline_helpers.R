make_group_folds <- function(data, v = 5, seed = 123) {
  set.seed(seed)
  patients <- unique(data$patient_id)
  fold_id <- sample(rep(seq_len(v), length.out = length(patients)))
  tibble::tibble(patient_id = patients, fold = fold_id)
}

split_by_fold <- function(data, fold_tbl, fold) {
  test_ids <- fold_tbl$patient_id[fold_tbl$fold == fold]
  list(
    train = dplyr::filter(data, !patient_id %in% test_ids),
    test  = dplyr::filter(data, patient_id %in% test_ids)
  )
}