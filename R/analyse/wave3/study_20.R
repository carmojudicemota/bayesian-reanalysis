study_20_claim_ids <- function() {
  c("study_20_claim_01", "study_20_claim_02")
}


study_20_normal_bf <- function(estimate, se, prior_sd = 0.5) {
  z <- estimate / se
  prior_ratio <- prior_sd^2 / se^2
  bf10_two_sided <- exp(-0.5 * log1p(prior_ratio) + 0.5 * z^2 * prior_ratio / (1 + prior_ratio))
  posterior_mean <- estimate * prior_sd^2 / (prior_sd^2 + se^2)
  posterior_sd <- sqrt(prior_sd^2 * se^2 / (prior_sd^2 + se^2))
  posterior_positive <- stats::pnorm(posterior_mean / posterior_sd)
  bf10_two_sided * 2 * posterior_positive
}


compute_study_20_bayes_factors <- function(claim, priors = NULL) {
  prior_sd <- 0.5
  bf10 <- study_20_normal_bf(
    estimate = as.numeric(claim$estimate),
    se = as.numeric(claim$se_estimate),
    prior_sd = prior_sd
  )
  wave3_row(
    claim = claim,
    bf10 = bf10,
    model_null = "contrast = 0",
    model_alt = "contrast > 0",
    bf_family = "normal_approximation",
    prior_family = "positive_half_normal",
    method = "wald_normal_approximation",
    prior_label = "primary"
  )
}
