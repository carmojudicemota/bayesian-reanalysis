library(RoBTT)
library(dplyr)
library(purrr)
library(haven)
library(tibble)

load_study_13_data <- function(outcome_col) {
  valid_outcomes <- c("subj_total", "obj_total")
  
  if (!outcome_col %in% valid_outcomes) {
    stop("Study 13 outcome must be 'subj_total' or 'obj_total'.", call. = FALSE)
  }
  
  raw <- haven::read_sav("data/raw/study_13/DATA_Cleaned_and_coded_for_condition.sav")
  
  data <- tibble(condition = as.numeric(raw$condition),outcome = as.numeric(raw[[outcome_col]])
  ) |>
    tidyr::drop_na(condition, outcome)
  
  groups <- list(intervention = data$outcome[data$condition == 1],control = data$outcome[data$condition == 0])
  
  if (length(groups$intervention) != 195L ||
      length(groups$control) != 303L) {
    stop(
      "Unexpected Study 13 sample sizes: intervention = ",
      length(groups$intervention),
      ", control = ",
      length(groups$control),
      ".",
      call. = FALSE
    )
  }
  
  groups
}

study_13_outcome <- function(claim_id) {
  switch(claim_id, 
         study_13_claim_01 = "subj_total",
         study_13_claim_02 = "obj_total",
         stop("Unsupported Study 13 claim: ", claim_id, call. = FALSE)
         )
}

study_13_prior_grid <- function(priors) {
  grid <- priors |>
    filter(.data$prior_family == "welch_averaged",
           .data$param == "delta_scale"
    ) |>
    transmute(prior_label = .data$prior_label,rscale = as.numeric(.data$value))
  
  if (
    nrow(grid) != 3L ||
    !setequal(grid$prior_label, c("narrow", "primary", "wide"))
  ) {
    stop("Study 13 requires narrow, primary and wide welch_averaged priors.",call. = FALSE)
  }
  
  grid
}

fit_study_13_model <- function(groups, rscale, seed = 2026) {
  RoBTT::RoBTT(
    x1 = groups$intervention,
    x2 = groups$control,
    prior_delta = RoBTT::prior("cauchy", list(0, rscale)),
    prior_rho = RoBTT::prior("beta", list(3, 3)),
    prior_nu = NULL,
    parallel = TRUE,
    seed = seed,
    control = RoBTT::set_control(adapt_delta = 0.95)
  )
}

extract_study_13_bf <- function(fit) {
  components <- summary(fit)$components
  if (!"Effect" %in% rownames(components)) {
    stop("RoBTT summary does not contain an Effect row.", call. = FALSE)
  }
  bf10 <- as.numeric(components["Effect", "inclusion_BF"])
  if (length(bf10) != 1L || !is.finite(bf10) || bf10 <= 0) {
    stop("RoBTT returned an invalid effect inclusion Bayes factor.", call. = FALSE)
  }
  
  bf10
}

run_study_13_claim <- function(claim_id, priors) {
  outcome_col <- study_13_outcome(claim_id)
  groups <- load_study_13_data(outcome_col)
  grid <- study_13_prior_grid(priors)
  
  purrr::map_dfr(seq_len(nrow(grid)), function(i) {
    prior_label <- grid$prior_label[[i]]
    rscale <- grid$rscale[[i]]
    
    message("Study 13: fitting ",claim_id," with ",prior_label," prior.")
    
    fit <- fit_study_13_model(groups = groups,rscale = rscale,seed = 2026 + i)
    fit_result <- extract_study_13_summary(fit)
    bf10 <- fit_result$bf10
    tibble(
      claim_id = claim_id,
      study_id = "study_13",
      prior_label = prior_label,
      rscale = rscale,
      bf10 = bf10,
      log_bf10 = log(bf10),
      log10_bf10 = log10(bf10),
      model_null = "delta = 0; equal- and unequal-variance normal models averaged",
      model_alt = paste0("delta ~ Cauchy(0, ",format(rscale, digits = 16),"); equal- and unequal-variance normal models averaged"),
      method = "RoBTT_variance_model_averaged",
      heterogeneity_bf = fit_result$heterogeneity_bf,
      delta_mean = fit_result$delta_mean,
      delta_median = fit_result$delta_median,
      delta_lower = fit_result$delta_lower,
      delta_upper = fit_result$delta_upper,
      rho_mean = fit_result$rho_mean
    )
  })
}

run_study_13_bayes_factors <- function(
    priors_path = "config/priors_wave2.csv") {
  priors <- readr::read_csv(priors_path, show_col_types = FALSE)
  
  bind_rows(
    run_study_13_claim("study_13_claim_01", priors),
    run_study_13_claim("study_13_claim_02", priors)
  )
}

extract_study_13_summary <- function(fit) {
  fit_summary <- summary(fit)
  components <- fit_summary$components
  estimates <- fit_summary$estimates
  
  tibble::tibble(
    bf10 = as.numeric(components["Effect", "inclusion_BF"]),
    heterogeneity_bf = as.numeric(components["Heterogeneity", "inclusion_BF"]),
    delta_mean = as.numeric(estimates["delta", "Mean"]),
    delta_median = as.numeric(estimates["delta", "Median"]),
    delta_lower = as.numeric(estimates["delta", "0.025"]),
    delta_upper = as.numeric(estimates["delta", "0.975"]),
    rho_mean = as.numeric(estimates["rho", "Mean"])
  )
}

validate_study_13_fit <- function(fit, claim_id, prior_label) {
  if (!all(fit$add_info$converged)) {
    stop("Study 13 model did not converge for ", claim_id, " under ", prior_label, ".", call. = FALSE)
  }
  
  invisible(TRUE)
}

compute_study_13_bayes_factors <- function(
    claim,
    priors,
    results_path = "outputs/intermediate/study_13_bayes_factors.csv") {
  
  if (!file.exists(results_path)) {
    stop("Study 13 cached results not found: ", results_path, ".", call. = FALSE)
  }
  
  cached <- readr::read_csv(results_path, show_col_types = FALSE) |>
    dplyr::filter(.data$claim_id == claim$claim_id) |>
    dplyr::mutate(
      prior_order = match(.data$prior_label, c("narrow", "primary", "wide"))
    ) |>
    dplyr::arrange(.data$prior_order)
  
  if (nrow(cached) != 3L || anyNA(cached$prior_order)) {
    stop("Study 13 requires three cached prior rows for ", claim$claim_id, ".", call. = FALSE)
  }
  
  purrr::pmap_dfr(cached, function(...) {
    row <- list(...)
    
    wave2_row(
      claim = claim,
      prior_label = row$prior_label,
      rscale = row$rscale,
      bf10 = row$bf10,
      log_bf10 = row$log_bf10,
      log10_bf10 = row$log10_bf10,
      model_null = row$model_null,
      model_alt = row$model_alt,
      bf_family = "welch_averaged",
      prior_family = "welch_averaged",
      method = row$method
    )
  })
}