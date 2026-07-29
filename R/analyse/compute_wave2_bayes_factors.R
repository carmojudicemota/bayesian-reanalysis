source("R/analyse/wave2/wave2_helpers.R")
source("R/analyse/wave2/factorial_anova.R")
source("R/analyse/wave2/mixed_anova.R")

source("R/analyse/wave2/study_29.R")
source("R/analyse/wave2/study_10.R")
source("R/analyse/wave2/study_60.R")
source("R/analyse/wave2/study_47.R")
source("R/analyse/wave2/study_35.R")
source("R/analyse/wave2/study_37_full_glmm.R")
source("R/analyse/wave2/study_06.R")
source("R/analyse/wave2/study_43.R")
source("R/analyse/wave2/study_13.R")
source("R/analyse/wave2/study_55.R")


compute_wave2_bayes_factors <- function(
    claims_path = "data/derived/claims.csv",
    priors_path = "config/priors_wave2.csv",
    output_path = "outputs/tables/bayes_factor_results_wave2.csv") {
  
  claims <- readr::read_csv(claims_path,show_col_types = FALSE)
  priors <- readr::read_csv(priors_path,show_col_types = FALSE )
  
  supported_claims <- c(
    "study_06_claim_01",
    "study_29_claim_02",
    "study_10_claim_01",
    "study_10_claim_02",
    "study_60_claim_01",
    "study_60_claim_02",
    "study_47_claim_02",
    "study_35_claim_01",
    "study_37_claim_01",
    "study_37_claim_02",
    "study_43_claim_01",
    "study_43_claim_02",
    "study_13_claim_01",
    "study_13_claim_02",
    "study_55_claim_01"
  )
  
  ready_claims <- claims |>
    dplyr::filter(.data$claim_id %in% supported_claims,
                  .data$status == "ready", .data$in_scope)
  
  results <- purrr::map_dfr(
    seq_len(nrow(ready_claims)),
    function(i) {
      claim <- as.list(ready_claims[i,,drop = FALSE])
      switch(
        claim$study_id,
        study_29 = compute_study_29_bayes_factors(claim = claim,priors = priors),
        study_10 = compute_study_10_bayes_factors(claim = claim,priors = priors),
        study_60 = compute_study_60_bayes_factors(claim = claim,priors = priors),
        study_47 = compute_study_47_bayes_factors(claim = claim,priors = priors),
        study_35 = compute_study_35_bayes_factors(claim = claim,priors = priors),
        study_37 = compute_study_37_bayes_factors(claim = claim,priors = priors),
        study_06 = compute_study_06_bayes_factors(claim = claim,priors = priors),
        study_43 = compute_study_43_bayes_factors(claim = claim, priors = priors),
        study_13 = compute_study_13_bayes_factors(claim = claim, priors = priors),
        study_55 = compute_study_55_bayes_factors(claim = claim, priors = priors),
        stop(
          "No Wave 2 implementation for ",
          claim$claim_id,
          ".",
          call. = FALSE
        )
      )
    }
  )
  
  if (nrow(ready_claims) == 0L) {results <- wave2_result_template() }
  dir.create(dirname(output_path),recursive = TRUE,showWarnings = FALSE)
  readr::write_csv(results,output_path,na = "")
  message(
    "Created Wave 2 output for ",
    dplyr::n_distinct(results$claim_id)," claims from ",
    dplyr::n_distinct(results$study_id)," studies."
    )
  invisible(results)
}