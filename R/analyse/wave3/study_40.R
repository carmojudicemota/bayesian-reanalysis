study_40_claim_ids <- function() {
  "study_40_claim_01"
}


load_study_40_data <- function(path = "data/raw/study_40/Wood_and_Cross_2024_Final_Data_deidentified.csv") {
  data <- utils::read.csv(path, stringsAsFactors = FALSE)
  data$opt_in <- factor(data$opt_in, levels = c(0, 1))
  outcomes <- paste0("p_tot", 1:5)
  data[outcomes] <- lapply(data[outcomes], as.numeric)
  scale_constant <- stats::sd(as.matrix(data[outcomes]), na.rm = TRUE)
  data[outcomes] <- data[outcomes] / scale_constant
  attr(data, "scale_constant") <- scale_constant
  data
}


study_40_model_7 <- function() {
  '
    i =~ 1*p_tot1 + 1*p_tot2 + 1*p_tot3 + 1*p_tot4 + 1*p_tot5
    s =~ 0*p_tot1 + m*p_tot2 + 1*p_tot3 + NA*p_tot4 + NA*p_tot5
  '
}


study_40_model_8 <- function() {
  '
    i =~ 1*p_tot1 + 1*p_tot2 + 1*p_tot3 + 1*p_tot4 + 1*p_tot5
    s =~ 0*p_tot1 + NA*p_tot2 + 1*p_tot3 + NA*p_tot4 + NA*p_tot5
  '
}


study_40_priors <- function(loading_scale = 0.5) {
  blavaan::dpriors(
    lambda = paste0("normal(1,", loading_scale, ")"),
    rho = "beta(2,2)"
  )
}


fit_study_40_model_7 <- function(data, burnin = 5000, sample = 10000, seed = 123,
                                 loading_scale = 0.5) {
  blavaan::bgrowth(
    model = study_40_model_7(),
    data = data,
    group = "opt_in",
    group.equal = c("residuals", "lv.variances", "lv.covariances"),
    cp = "srs",
    dp = study_40_priors(loading_scale),
    n.chains = 4,
    burnin = burnin,
    sample = sample,
    adapt = 2000,
    seed = seed,
    bcontrol = list(control = list(adapt_delta = 0.999, max_treedepth = 15))
  )
}


fit_study_40_model_8 <- function(data, burnin = 5000, sample = 10000, seed = 124,
                                 loading_scale = 0.5) {
  blavaan::bgrowth(
    model = study_40_model_8(),
    data = data,
    group = "opt_in",
    group.equal = c("loadings", "residuals", "lv.variances", "lv.covariances"),
    cp = "srs",
    dp = study_40_priors(loading_scale),
    n.chains = 4,
    burnin = burnin,
    sample = sample,
    adapt = 2000,
    seed = seed,
    bcontrol = list(control = list(adapt_delta = 0.999, max_treedepth = 15))
  )
}


study_40_diagnostics <- function(fit) {
  rhat <- blavaan::blavInspect(fit, "rhat")
  ess <- blavaan::blavInspect(fit, "neff")
  list(
    rhat_max = max(rhat, na.rm = TRUE),
    ess_bulk_min = min(ess, na.rm = TRUE),
    converged = max(rhat, na.rm = TRUE) < 1.01 && min(ess, na.rm = TRUE) > 400
  )
}


study_40_cache_path <- function() {
  "outputs/intermediate/study_40_bayes_factors.csv"
}


study_40_write_cache <- function(fit7, fit8, loading_scale = 0.5, cache_path = study_40_cache_path()) {
  diag7 <- study_40_diagnostics(fit7)
  diag8 <- study_40_diagnostics(fit8)
  log_bf_7_8 <- as.numeric(lavaan::fitMeasures(fit7, "margloglik") -
                             lavaan::fitMeasures(fit8, "margloglik"))
  result <- data.frame(
    claim_id = "study_40_claim_01",
    bf10 = exp(log_bf_7_8),
    log_bf10 = log_bf_7_8,
    rhat_max_m7 = diag7$rhat_max,
    rhat_max_m8 = diag8$rhat_max,
    loading_scale = loading_scale,
    method = "blavaan_bgrowth_margloglik",
    generated_at = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result, cache_path, row.names = FALSE)
  invisible(result)
}


run_study_40_bayes_factors <- function(cache_path = study_40_cache_path(), loading_scale = 0.5, seed = 123) {
  data <- load_study_40_data()
  fit7 <- fit_study_40_model_7(data, seed = seed, loading_scale = loading_scale)
  fit8 <- fit_study_40_model_8(data, seed = seed + 1, loading_scale = loading_scale)
  diag7 <- study_40_diagnostics(fit7)
  diag8 <- study_40_diagnostics(fit8)
  if (!diag7$converged || !diag8$converged) {
    stop("Study 40 fits did not converge; not caching. rhat_max M7=", round(diag7$rhat_max, 4),
         " M8=", round(diag8$rhat_max, 4), call. = FALSE)
  }
  saveRDS(fit7, file.path(dirname(cache_path), "study_40_fit7.rds"))
  saveRDS(fit8, file.path(dirname(cache_path), "study_40_fit8.rds"))
  study_40_write_cache(fit7, fit8, loading_scale, cache_path)
}


study_40_cache_from_saved <- function(fit7_path = "outputs/intermediate/study_40_fit7.rds",
                                      fit8_path = "outputs/intermediate/study_40_fit8.rds",
                                      loading_scale = 0.5, cache_path = study_40_cache_path()) {
  study_40_write_cache(readRDS(fit7_path), readRDS(fit8_path), loading_scale, cache_path)
}


compute_study_40_bayes_factors <- function(claim, priors = NULL, cache_path = study_40_cache_path()) {
  if (!file.exists(cache_path)) {
    stop("Study 40 cache not found at ", cache_path, ". Run run_study_40_bayes_factors() first.", call. = FALSE)
  }
  cached <- utils::read.csv(cache_path, stringsAsFactors = FALSE)
  row <- cached[cached$claim_id == claim$claim_id, ]
  if (nrow(row) != 1) {
    stop("No cached Study 40 result for ", claim$claim_id, ".", call. = FALSE)
  }
  wave3_row(
    claim = claim,
    bf10 = row$bf10,
    model_null = "Model 8: loadings constrained equal across groups",
    model_alt = "Model 7: group-specific post-intervention loadings",
    bf_family = "sem_marginal_likelihood",
    prior_family = "blavaan_default",
    method = "blavaan_bgrowth_blavCompare",
    prior_label = "primary"
  )
}
