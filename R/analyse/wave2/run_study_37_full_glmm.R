study_37_focal_coefficient <- function(claim_id) {
  switch(
    claim_id,
    study_37_claim_01 = "HUMAN_family1",
    study_37_claim_02 = "synchronous11",
    stop("Unknown Study 37 claim: ", claim_id, call. = FALSE)
  )
}

study_37_prior_values <- function(focal_sd) {
  list(nuisance_sd = 1, intercept_sd = 2.5, random_sd = 1, focal_sd = focal_sd)
}

study_37_prior_grid <- function() {
  tibble::tibble(
    prior_label = c("narrow", "primary", "wide"),
    focal_sd = c(0.3, 0.6, 1)
  )
}

fit_study_37_model <- function(formula, data, priors, seed,
                               iter = 6000, warmup = 2000, chains = 4, cores = 4) {
  brms::brm(
    formula = formula,
    data = data,
    family = brms::bernoulli(),
    prior = priors,
    iter = iter,
    warmup = warmup,
    chains = chains,
    cores = cores,
    seed = seed,
    save_pars = brms::save_pars(all = TRUE),
    refresh = 0
  )
}

run_study_37_full_glmm <- function(output_path = "outputs/tables/study_37_full_glmm_bayes_factors.csv",
                                   repetitions = 5L, cores = 4L, base_seed = 123L) {
  data <- load_study_37_full_glmm_data()
  claims <- c("study_37_claim_01", "study_37_claim_02")
  grid <- study_37_prior_grid()
  rows <- purrr::map_dfr(claims, function(claim_id) {
    focal_coef <- study_37_focal_coefficient(claim_id)
    full_formula <- study_37_full_formula()
    null_formula <- study_37_null_formula(claim_id)
    purrr::map_dfr(seq_len(nrow(grid)), function(i) {
      prior_label <- grid$prior_label[[i]]
      focal_sd <- grid$focal_sd[[i]]
      prior_values <- study_37_prior_values(focal_sd)
      seed_full <- base_seed + 2L * i
      seed_null <- base_seed + 2L * i + 1L
      full_fit <- fit_study_37_model(full_formula, data,
                                     study_37_model_priors(claim_id, prior_values, TRUE),
                                     seed_full, cores = cores)
      null_fit <- fit_study_37_model(null_formula, data,
                                     study_37_model_priors(claim_id, prior_values, FALSE),
                                     seed_null, cores = cores)
      set.seed(seed_full)
      bridge <- bridge_study_37_pair(list(full = full_fit, null = null_fit),
                                     repetitions = repetitions, cores = cores)
      central <- stats::median(bridge$log10_bf10)
      tibble::tibble(
        claim_id = claim_id,
        study_id = "study_37",
        prior_label = prior_label,
        rscale = focal_sd,
        bf10 = 10^central,
        log_bf10 = central * log(10),
        log10_bf10 = central,
        bf_error = NA_real_,
        bridge_sd_log10 = stats::sd(bridge$log10_bf10),
        bridge_span_log10 = max(bridge$log10_bf10) - min(bridge$log10_bf10),
        bridge_repetitions = nrow(bridge),
        model_null = sprintf("full Bayesian GLMM without %s", focal_coef),
        model_alt = sprintf("full Bayesian GLMM with %s", focal_coef),
        method = "full_bayesian_glmm_bridge_sampling"
      )
    })
  })
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(rows, output_path, na = "")
  message(sprintf("Wrote %d Study 37 rows to %s", nrow(rows), output_path))
  invisible(rows)
}
