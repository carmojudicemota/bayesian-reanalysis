source("R/analyse/wave3/wave3_helpers.R")

source("R/analyse/wave3/study_03.R")


compute_wave3_bayes_factors <- function(claims_path = "data/derived/claims.csv",output_path = "outputs/tables/bayes_factor_results_wave3.csv",iter = 100000) {
  
  claims <- readr::read_csv(claims_path,show_col_types = FALSE)
  supported_claims <- study_03_claim_ids()
  ready_claims <- claims |>
    dplyr::filter(.data$claim_id %in% supported_claims,.data$status == "ready",.data$in_scope)

  results <- purrr::map_dfr(
    seq_len(nrow(ready_claims)),
    function(i) {
      claim <- as.list(ready_claims[i, , drop = FALSE])
      switch(
        claim$study_id,
        study_03 = compute_study_03_bayes_factors(claim = claim,iter = iter),
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