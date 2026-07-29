wave3_result_template <- function() {
  tibble::tibble(
    claim_id = character(),
    study_id = character(),
    stat_test = character(),
    bf_family = character(),
    design = character(),
    bf_sidedness = character(),
    bf_direction = character(),
    observed_sign = character(),
    direction_matches_observed = logical(),
    prior_label = character(),
    rscale = double(),
    t_for_bf = double(),
    df_for_bf = double(),
    r_value = double(),
    n1 = double(),
    n2 = double(),
    n_total = double(),
    p_value = double(),
    bf10 = double(),
    log_bf10 = double(),
    log10_bf10 = double(),
    bf_error = double(),
    prior_family = character(),
    model_null = character(),
    model_alt = character(),
    method = character()
  )
}


wave3_row <- function(
    claim,
    bf10,
    model_null,
    model_alt,
    bf_family,
    prior_family,
    method,
    prior_label = "primary",
    bf_error = NA_real_) {
  
  bf10 <- as.numeric(bf10)
  
  tibble::tibble(
    claim_id = claim$claim_id,
    study_id = claim$study_id,
    stat_test = claim$frequentist_test,
    bf_family = bf_family,
    design = "wave3",
    bf_sidedness = "two_sided",
    bf_direction = NA_character_,
    observed_sign = NA_character_,
    direction_matches_observed = NA,
    prior_label = prior_label,
    rscale = NA_real_,
    t_for_bf = NA_real_,
    df_for_bf = NA_real_,
    r_value = NA_real_,
    n1 = as.numeric(claim$n1),
    n2 = as.numeric(claim$n2),
    n_total = as.numeric(claim$n_total),
    p_value = as.numeric(claim$p_value),
    bf10 = bf10,
    log_bf10 = log(bf10),
    log10_bf10 = log10(bf10),
    bf_error = as.numeric(bf_error),
    prior_family = prior_family,
    model_null = model_null,
    model_alt = model_alt,
    method = method
  )
}


wave3_manova_diagnostics <- function(model) {
  resid <- stats::residuals(model)
  shapiro_p <- apply(resid, 2, function(z) stats::shapiro.test(z)$p.value)
  n <- nrow(resid)
  list(
    n = n,
    n_outcomes = ncol(resid),
    n_location_parameters = nrow(stats::coef(model)),
    fraction_b = nrow(stats::coef(model)) / n,
    min_shapiro_p = min(shapiro_p),
    shapiro_p = shapiro_p
  )
}


wave3_check_assumptions <- function(model, claim_id, min_shapiro = 0.001) {
  d <- wave3_manova_diagnostics(model)
  message(sprintf(
    "%s: N=%d, outcomes=%d, fraction b~%.3f, min Shapiro p=%.4g",
    claim_id, d$n, d$n_outcomes, d$fraction_b, d$min_shapiro_p
  ))
  if (d$min_shapiro_p < min_shapiro) {
    warning(sprintf(
      "%s: residual normality is doubtful (min Shapiro p=%.4g); consider the Student-t robustness fit.",
      claim_id, d$min_shapiro_p
    ), call. = FALSE)
  }
  invisible(d)
}


wave3_manova_diagnostics <- function(model) {
  resid <- stats::residuals(model)
  shapiro_p <- apply(resid, 2, function(z) stats::shapiro.test(z)$p.value)
  n <- nrow(resid)
  box_m_p <- NA_real_
  list(
    n = n,
    n_outcomes = ncol(resid),
    n_location_parameters = nrow(stats::coef(model)),
    fraction_b = nrow(stats::coef(model)) / n,
    min_shapiro_p = min(shapiro_p),
    shapiro_p = shapiro_p,
    box_m_p = box_m_p
  )
}


wave3_check_assumptions <- function(model, claim_id, min_shapiro = 0.001) {
  d <- wave3_manova_diagnostics(model)
  message(sprintf(
    "%s: N=%d, outcomes=%d, fraction b~%.3f, min Shapiro p=%.4g",
    claim_id, d$n, d$n_outcomes, d$fraction_b, d$min_shapiro_p
  ))
  if (d$min_shapiro_p < min_shapiro) {
    warning(sprintf(
      "%s: residual normality is doubtful (min Shapiro p=%.4g); consider the Student-t robustness fit.",
      claim_id, d$min_shapiro_p
    ), call. = FALSE)
  }
  invisible(d)
}