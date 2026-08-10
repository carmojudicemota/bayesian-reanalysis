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


study_40_loading_scales <- function() {
  c(narrow = 0.25, primary = 0.5, wide = 1.0)
}


study_40_cache_row <- function(fit7, fit8, loading_scale, prior_label) {
  diag7 <- study_40_diagnostics(fit7)
  diag8 <- study_40_diagnostics(fit8)
  log_bf_7_8 <- as.numeric(lavaan::fitMeasures(fit7, "margloglik") -
                             lavaan::fitMeasures(fit8, "margloglik"))
  data.frame(
    claim_id = "study_40_claim_01",
    prior_label = prior_label,
    loading_scale = loading_scale,
    bf10 = exp(log_bf_7_8),
    log_bf10 = log_bf_7_8,
    rhat_max_m7 = diag7$rhat_max,
    rhat_max_m8 = diag8$rhat_max,
    method = "blavaan_bgrowth_margloglik",
    generated_at = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
}


study_40_fit_scale <- function(loading_scale, prior_label, seed = 123,
                               cache_path = study_40_cache_path(), save_fits = TRUE) {
  data <- load_study_40_data()
  fit7 <- fit_study_40_model_7(data, seed = seed, loading_scale = loading_scale)
  fit8 <- fit_study_40_model_8(data, seed = seed + 1, loading_scale = loading_scale)
  diag7 <- study_40_diagnostics(fit7)
  diag8 <- study_40_diagnostics(fit8)
  if (!diag7$converged || !diag8$converged) {
    stop("Study 40 fits did not converge at loading_scale ", loading_scale,
         "; rhat_max M7=", round(diag7$rhat_max, 4), " M8=", round(diag8$rhat_max, 4), call. = FALSE)
  }
  if (save_fits) {
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(fit7, file.path(dirname(cache_path), paste0("study_40_fit7_", prior_label, ".rds")))
    saveRDS(fit8, file.path(dirname(cache_path), paste0("study_40_fit8_", prior_label, ".rds")))
  }
  study_40_cache_row(fit7, fit8, loading_scale, prior_label)
}


run_study_40_bayes_factors <- function(cache_path = study_40_cache_path(),
                                       loading_scales = study_40_loading_scales(), seed = 123) {
  rows <- purrr::map_dfr(seq_along(loading_scales), function(i) {
    study_40_fit_scale(as.numeric(loading_scales[i]), names(loading_scales)[i],
                       seed = seed, cache_path = cache_path)
  })
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(rows, cache_path, row.names = FALSE)
  invisible(rows)
}


study_40_cache_from_saved <- function(loading_scales = study_40_loading_scales(),
                                      cache_path = study_40_cache_path()) {
  rows <- purrr::map_dfr(seq_along(loading_scales), function(i) {
    prior_label <- names(loading_scales)[i]
    fit7 <- readRDS(file.path(dirname(cache_path), paste0("study_40_fit7_", prior_label, ".rds")))
    fit8 <- readRDS(file.path(dirname(cache_path), paste0("study_40_fit8_", prior_label, ".rds")))
    study_40_cache_row(fit7, fit8, as.numeric(loading_scales[i]), prior_label)
  })
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(rows, cache_path, row.names = FALSE)
  invisible(rows)
}


study_40_predictive_comparison <- function(
    fit7_path = "outputs/intermediate/study_40_fit7_primary.rds",
    fit8_path = "outputs/intermediate/study_40_fit8_primary.rds",
    output_path = "outputs/intermediate/study_40_predictive_comparison.csv") {
  if (!file.exists(fit7_path) || !file.exists(fit8_path)) {
    message("Study 40 primary fits not found; running the loading sweep to produce them.")
    run_study_40_bayes_factors()
  }
  fit7 <- readRDS(fit7_path)
  fit8 <- readRDS(fit8_path)
  measures7 <- lavaan::fitMeasures(fit7, c("waic", "looic"))
  measures8 <- lavaan::fitMeasures(fit8, c("waic", "looic"))
  waic7 <- as.numeric(measures7[["waic"]])
  waic8 <- as.numeric(measures8[["waic"]])
  looic7 <- as.numeric(measures7[["looic"]])
  looic8 <- as.numeric(measures8[["looic"]])
  row <- tibble::tibble(
    claim_id = "study_40_claim_01",
    model_alt = "Model 7: group-specific post-intervention loadings",
    model_null = "Model 8: loadings constrained equal across groups",
    waic_alt = waic7,
    waic_null = waic8,
    waic_diff_alt_minus_null = waic7 - waic8,
    looic_alt = looic7,
    looic_null = looic8,
    looic_diff_alt_minus_null = looic7 - looic8,
    elpd_diff_alt_minus_null = -(looic7 - looic8) / 2
  )
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(row, output_path, na = "")
  print(as.data.frame(row))
  invisible(row)
}

compute_study_40_bayes_factors <- function(claim, priors = NULL, cache_path = study_40_cache_path()) {
  if (!file.exists(cache_path)) {
    stop("Study 40 cache not found at ", cache_path, ". Run run_study_40_bayes_factors() first.", call. = FALSE)
  }
  cached <- utils::read.csv(cache_path, stringsAsFactors = FALSE)
  rows <- cached[cached$claim_id == claim$claim_id, , drop = FALSE]
  if (nrow(rows) < 1) {
    stop("No cached Study 40 result for ", claim$claim_id, ".", call. = FALSE)
  }
  if (is.null(rows$prior_label)) {
    rows$prior_label <- "primary"
  }
  purrr::map_dfr(seq_len(nrow(rows)), function(i) {
    wave3_row(
      claim = claim,
      bf10 = rows$bf10[i],
      model_null = "Model 8: loadings constrained equal across groups",
      model_alt = "Model 7: group-specific post-intervention loadings",
      bf_family = "sem_marginal_likelihood",
      prior_family = "blavaan_default",
      method = "blavaan_bgrowth_margloglik",
      prior_label = rows$prior_label[i]
    )
  })
}
