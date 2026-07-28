study_03_claim_ids <- function() {c("study_03_claim_01", "study_03_claim_02")}

load_study_03_wave3_data <- function(
    path = "data/raw/study_03/Snapshot_MASTER_noID_cleaned_03.11.24.sav") {
  
  raw <- haven::read_sav(path)
  clean_numeric <- function(x) {
    x <- haven::zap_labels(x)
    if (is.factor(x)) x <- as.character(x)
    suppressWarnings(as.numeric(x))
  }
  
  tibble::tibble(
    TBC_COMP = clean_numeric(raw$TBC_COMP),
    TBC_CARE = clean_numeric(raw$TBC_CARE),
    percept = clean_numeric(raw$percept),
    snap_syll = clean_numeric(raw$snap_syll),
    snap_inst = clean_numeric(raw$snap_inst)
  ) |>
    dplyr::filter(
      stats::complete.cases(
        TBC_COMP,
        TBC_CARE,
        percept,
        snap_syll,
        snap_inst
      )
    ) |>
    dplyr::mutate(
      syllabus_snapshot = ifelse(snap_syll == 1, 1, -1),
      instructor_snapshot = ifelse(snap_inst == 1, 1, -1),
      syllabus_x_instructor = syllabus_snapshot * instructor_snapshot
    )
}

fit_study_03_multivariate_model <- function(data) {
  stats::lm(cbind(TBC_COMP, TBC_CARE, percept) ~ syllabus_snapshot +
              instructor_snapshot + syllabus_x_instructor,
            data = data
  )
}



