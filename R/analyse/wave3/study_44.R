study_44_claim_ids <- function() {
  c("study_44_claim_01","study_44_claim_02")
}


load_study_44_data <- function(path = "data/raw/study_44/Untitled2.sav") {
  raw <- haven::read_sav(path)
  
  as_number <- function(x) {suppressWarnings(as.numeric(as.character(x)))}
  
  data <- raw |>
    dplyr::transmute(
      experiment = as_number(METADATA),
      participant_status_code = as_number(Group),
      condition_code = as_number(V6),
      effective = as_number(OUTCOMES1),
      well_liked = as_number(V8),
      clear_expectations_reversed = as_number(V10),
      helpful = as_number(V11),
      challenging_content = as_number(V12),
      learning_environment = as_number(V13),
      prepared_reversed = as_number(V15)
    ) |>
    dplyr::filter(experiment == 1) |>
    dplyr::mutate(
      participant_status = factor(
        participant_status_code,
        levels = c(1, 2),
        labels = c("student","faculty")
      ),
      rating_magnitude = factor(
        dplyr::case_when(
          condition_code %in% c(1, 3) ~ "high",
          condition_code %in% c(2, 4) ~ "low",
          TRUE ~ NA_character_
        ),
        levels = c("high","low")
      ),
      quizzing_comment = factor(
        dplyr::case_when(
          condition_code %in% c(1, 2) ~ "quizzes",
          condition_code %in% c(3, 4) ~ "no_quizzes",
          TRUE ~ NA_character_
        ),
        levels = c("quizzes","no_quizzes")
      )
    ) |>
    dplyr::select(
      participant_status,
      rating_magnitude,
      quizzing_comment,
      effective,
      well_liked,
      clear_expectations_reversed,
      helpful,
      challenging_content,
      learning_environment,
      prepared_reversed
    ) |>
    tidyr::drop_na()
  
  contrasts(data$participant_status) <- stats::contr.sum(2)
  contrasts(data$rating_magnitude) <- stats::contr.sum(2)
  contrasts(data$quizzing_comment) <- stats::contr.sum(2)
  
  data
}


fit_study_44_model <- function(data) {
  stats::lm(
    cbind(
      effective,
      well_liked,
      clear_expectations_reversed,
      helpful,
      challenging_content,
      learning_environment,
      prepared_reversed
    ) ~ participant_status * rating_magnitude * quizzing_comment,
    data = data
  )
}

study_44_outcomes <- function() {
  c(
    "effective",
    "well_liked",
    "clear_expectations_reversed",
    "helpful",
    "challenging_content",
    "learning_environment",
    "prepared_reversed"
  )
}


study_44_hypothesis <- function(claim_id) {
  coefficient <- switch(
    claim_id,
    study_44_claim_01 = "rating_magnitude1",
    study_44_claim_02 = "quizzing_comment1",
    stop("Unknown Study 44 claim: ",claim_id,call. = FALSE)
  )
  paste0(coefficient,"_on_",study_44_outcomes()," = 0",collapse = " & ")
}


study_44_model_labels <- function(claim_id) {
  switch(
    claim_id,
    study_44_claim_01 = list(
      null = "all rating-magnitude coefficients = 0",
      alternative = "at least one rating-magnitude coefficient != 0"
    ),
    study_44_claim_02 = list(
      null = "all quizzing-comment coefficients = 0",
      alternative = "at least one quizzing-comment coefficient != 0"
    ),
    stop("Unknown Study 44 claim: ",claim_id,call. = FALSE)
  )
}

compute_study_44_bayes_factors <- function(claim,priors = NULL,iter = 100000) {
  data <- load_study_44_data()
  model <- fit_study_44_model(data)
  hypothesis <- study_44_hypothesis(claim$claim_id)
  
  set.seed(123)
  
  result <- BFpack::BF(model,hypothesis = hypothesis,iter = iter)
  
  bf10 <- as.numeric(result$BFmatrix_confirmatory["H2", "H1"])
  
  labels <- study_44_model_labels(claim$claim_id)
  
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


