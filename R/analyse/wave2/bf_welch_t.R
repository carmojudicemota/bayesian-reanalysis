library(RoBTT)
library(dplyr)
library(purrr)

run_family_C <- function(claim, priors) {
  gp <- load_wave2_data(claim$study_id, outcome_col = claim$outcome_column)
  
  grid <- priors |>
    filter(prior_family == "welch_averaged", param == "delta_scale")
  
  pmap_dfr(grid, function(prior_label, value, ...) {
    fit <- RoBTT::RoBTT(
      x1 = gp$x1, x2 = gp$x2,
      prior_delta = RoBTT::prior("cauchy", list(location = 0, scale = value)),
      prior_rho = RoBTT::prior("beta", list(alpha = 1, beta = 1)),
      chains = 4, iter = 8000, warmup = 2000, seed = 2026
    )
    
    fit_summary <- summary(fit)
    diagnostics <- fit_summary$diagnostics
    
    record_mcmc_diagnostics(
      claim$claim_id,
      diagnostics$rhat_max,
      diagnostics$ess_min,
      diagnostics$divergences,
      diagnostics$ppc_ok
    )
    
    wave2_row(
      claim = claim,
      prior_label = prior_label,
      rscale = value,
      bf10 = fit_summary$estimates["delta", "BF"],
      bf_family = "welch_averaged",
      model_null = "delta = 0; equal- and unequal-variance models averaged",
      model_alt = "delta ~ Cauchy(0, r); equal- and unequal-variance models averaged"
    )
  })
}