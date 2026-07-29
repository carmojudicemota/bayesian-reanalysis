study_03_claim_ids <- function() {
  c("study_03_claim_01", "study_03_claim_02")
}


load_study_03_data <- function(path = NULL) {
  candidates <- c(
    "data/raw/study_03/Snapshot_MASTER_noID_cleaned_03.11.24.sav",
    "data/raw/study_03/Snapshot MASTER_noID cleaned 03.11.24.sav"
  )
  if (is.null(path)) {
    path <- candidates[file.exists(candidates)][1]
  }
  if (is.na(path) || !file.exists(path)) {
    stop("Study 03 data file not found.", call. = FALSE)
  }
  raw <- haven::read_sav(path)
  as_number <- function(x) {
    x <- haven::zap_labels(x)
    if (is.factor(x)) x <- as.character(x)
    suppressWarnings(as.numeric(x))
  }
  raw |>
    dplyr::transmute(
      teacher_competence = as_number(TBC_COMP),
      teacher_care = as_number(TBC_CARE),
      general_perception = as_number(percept),
      syllabus_snapshot = ifelse(as_number(snap_syll) == 1, 1, -1),
      instructor_snapshot = ifelse(as_number(snap_inst) == 1, 1, -1)
    ) |>
    tidyr::drop_na()
}


study_03_outcomes <- function() {
  c("teacher_competence", "teacher_care", "general_perception")
}


fit_study_03_model <- function(data) {
  stats::lm(
    cbind(teacher_competence, teacher_care, general_perception) ~
      syllabus_snapshot * instructor_snapshot,
    data = data
  )
}


study_03_hypothesis <- function(claim_id) {
  term <- switch(
    claim_id,
    study_03_claim_01 = "syllabus_snapshot:instructor_snapshot",
    study_03_claim_02 = "syllabus_snapshot",
    stop("Unknown Study 03 claim: ", claim_id, call. = FALSE)
  )
  paste0(term, "_on_", study_03_outcomes(), " = 0", collapse = " & ")
}


study_03_model_labels <- function(claim_id) {
  switch(
    claim_id,
    study_03_claim_01 = list(
      null = "syllabus-by-instructor interaction block = 0",
      alternative = "at least one interaction coefficient != 0"
    ),
    study_03_claim_02 = list(
      null = "all syllabus-snapshot coefficients = 0",
      alternative = "at least one syllabus-snapshot coefficient != 0"
    ),
    stop("Unknown Study 03 claim: ", claim_id, call. = FALSE)
  )
}


compute_study_03_bayes_factors <- function(claim, priors = NULL, iter = 100000) {
  data <- load_study_03_data()
  model <- fit_study_03_model(data)
  wave3_check_assumptions(model, claim$claim_id)
  set.seed(123)
  result <- BFpack::BF(model, hypothesis = study_03_hypothesis(claim$claim_id), iter = iter)
  bf10 <- as.numeric(result$BFmatrix_confirmatory["H2", "H1"])
  labels <- study_03_model_labels(claim$claim_id)
  wave3_row(
    claim = claim,
    bf10 = bf10,
    model_null = labels$null,
    model_alt = labels$alternative,
    bf_family = "generalized_fractional",
    prior_family = "fractional",
    method = "BFpack_generalized_fractional",
    prior_label = "primary"
  )
}
