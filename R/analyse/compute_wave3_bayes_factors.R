source("R/analyse/wave3/wave3_helpers.R")

source("R/analyse/wave3/study_03.R")
source("R/analyse/wave3/study_18.R")
source("R/analyse/wave3/study_44.R")
source("R/analyse/wave3/study_20.R")
source("R/analyse/wave3/study_40.R")


compute_wave3_bayes_factors <- function(claims_path = "data/derived/claims.csv",output_path = "outputs/tables/bayes_factor_results_wave3.csv",iter = 100000) {
  
  claims <- readr::read_csv(claims_path,show_col_types = FALSE)
  supported_claims <- c(
    study_03_claim_ids(),
    study_18_claim_ids(),
    study_44_claim_ids(),
    study_20_claim_ids(),
    study_40_claim_ids(),
    study_39_claim_ids()
  )
  ready_claims <- claims |>
    dplyr::filter(.data$claim_id %in% supported_claims,
                  .data$status == "ready")

  results <- purrr::map_dfr(
    seq_len(nrow(ready_claims)),
    function(i) {
      claim <- as.list(ready_claims[i, , drop = FALSE])
      switch(claim$study_id,
             study_03 = compute_study_03_bayes_factors(claim = claim),
             study_18 = compute_study_18_bayes_factors(claim = claim),
             study_44 = compute_study_44_bayes_factors(claim = claim),
             study_20 = compute_study_20_bayes_factors(claim = claim),
             study_40 = compute_study_40_bayes_factors(claim = claim),
             study_39 = compute_study_39_bayes_factors(claim = claim, n_samples = 20000, n_burnin = 5000),
             stop("No Wave 3 implementation for ",claim$claim_id,".",call. = FALSE)
             )
    }
  )
  
  if (nrow(ready_claims) == 0L) {results <- wave3_result_template()}
  dir.create(dirname(output_path),recursive = TRUE,showWarnings = FALSE)
  readr::write_csv(results,output_path, na = "")
  message("Created Wave 3 output for ", dplyr::n_distinct(results$claim_id),
          " claims from ", dplyr::n_distinct(results$study_id)," studies.")
  
  invisible(results)
}