study_03_claim_ids <- function() {
  c(
    "study_03_claim_01",
    "study_03_claim_02"
  )
}


load_study_03_data <- function(
    path = "data/raw/study_03/Snapshot_MASTER_noID_cleaned_03.11.24.sav") {
  
  raw <- haven::read_sav(path)
  
  to_numeric <- function(x) {
    as.numeric(haven::zap_labels(x))
  }
  
  tibble::tibble(
    TBC_COMP = to_numeric(raw$TBC_COMP),
    TBC_CARE = to_numeric(raw$TBC_CARE),
    percept = to_numeric(raw$percept),
    snap_syll = to_numeric(raw$snap_syll),
    snap_inst = to_numeric(raw$snap_inst)
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
      syllabus = ifelse(snap_syll == 1, 1, -1),
      instructor = ifelse(snap_inst == 1, 1, -1),
      syllabus_instructor = syllabus * instructor
    )
}


fit_study_03_model <- function(data) {
  stats::lm(
    cbind(TBC_COMP, TBC_CARE, percept) ~
      syllabus +
      instructor +
      syllabus_instructor,
    data = data
  )
}


study_03_hypothesis <- function(claim_id) {
  switch(
    claim_id,
    
    study_03_claim_01 = paste(
      "syllabus_instructor_on_TBC_COMP = 0",
      "syllabus_instructor_on_TBC_CARE = 0",
      "syllabus_instructor_on_percept = 0",
      sep = " & "
    ),
    
    study_03_claim_02 = paste(
      "syllabus_on_TBC_COMP = 0",
      "syllabus_on_TBC_CARE = 0",
      "syllabus_on_percept = 0",
      sep = " & "
    ),
    
    stop(
      "Unknown Study 03 claim: ",
      claim_id,
      call. = FALSE
    )
  )
}


study_03_model_labels <- function(claim_id) {
  switch(
    claim_id,
    
    study_03_claim_01 = list(
      null = "joint interaction coefficients = 0",
      alternative = "at least one interaction coefficient != 0"
    ),
    
    study_03_claim_02 = list(
      null = "joint syllabus coefficients = 0",
      alternative = "at least one syllabus coefficient != 0"
    ),
    
    stop(
      "Unknown Study 03 claim: ",
      claim_id,
      call. = FALSE
    )
  )
}


compute_study_03_bayes_factors <- function(
    claim,
    priors = NULL,
    iter = 100000) {
  
  data <- load_study_03_data()
  model <- fit_study_03_model(data)
  
  hypothesis <- study_03_hypothesis(
    claim$claim_id
  )
  
  set.seed(123)
  
  result <- BFpack::BF(
    model,
    hypothesis = hypothesis,
    iter = iter
  )
  
  bf01 <- as.numeric(
    result$BFmatrix_confirmatory["H1", "H2"]
  )
  
  labels <- study_03_model_labels(
    claim$claim_id
  )
  
  wave3_row(
    claim = claim,
    bf10 = 1 / bf01,
    model_null = labels$null,
    model_alt = labels$alternative,
    bf_family = "generalized_fractional",
    prior_family = "fractional",
    method = "BFpack_generalized_fractional",
    prior_label = "primary"
  )
}