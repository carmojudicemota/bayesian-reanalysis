library(brms)
library(bridgesampling)
library(dplyr)
library(purrr)
library(readr)
library(tibble)
library(tidyr)

load_study_55_data <- function(
    path = "data/raw/study_55/MQ_Exp_1_data_removed_failed_mastery_Ss_and_timed_out.csv") {
  
  raw <- readr::read_csv(path, show_col_types = FALSE)
  data <- raw |>
    dplyr::transmute(
      condition = trimws(as.character(.data$Condition)),
      outcome = as.numeric(.data[["Final Test Score"]])
    ) |>
    tidyr::drop_na(condition, outcome) |>
    dplyr::filter(.data$condition %in% c("one", "two", "three", "mastery")) |>
    dplyr::mutate(
      condition = factor(.data$condition, levels = c("one", "two", "three", "mastery"))
    )
  
  if (nrow(data) != 154L) {
    stop("Unexpected Study 55 sample size: ", nrow(data), ".", call. = FALSE)
  }
  
  outcome_sd <- stats::sd(data$outcome)
  if (!is.finite(outcome_sd) || outcome_sd <= 0) {
    stop("Study 55 outcome standard deviation is invalid.", call. = FALSE)
  }
  
  contrast_matrix <- stats::contr.helmert(4)
  contrast_matrix <- sweep(contrast_matrix, 2, sqrt(colSums(contrast_matrix^2)), "/")
  design <- contrast_matrix[as.integer(data$condition), , drop = FALSE]
  colnames(design) <- c("contrast_1", "contrast_2", "contrast_3")
  dplyr::bind_cols(
    data |>
      dplyr::mutate(
        outcome_z = (.data$outcome - mean(.data$outcome)) / outcome_sd,
        Intercept = 1
        ), 
    tibble::as_tibble(design)
    )
}

validate_study_55_reconstruction <- function(data) {
  welch <- stats::oneway.test(outcome ~ condition, data = data, var.equal = FALSE)
  
  validation <- tibble::tibble(
    f_value = unname(welch$statistic),
    df1 = unname(welch$parameter[[1]]),
    df2 = unname(welch$parameter[[2]]),
    p_value = welch$p.value,
    n_total = nrow(data)
  )
  
  if (
    abs(validation$f_value - 35.64) > 0.1 ||
    abs(validation$df1 - 3) > 0.01 ||
    abs(validation$df2 - 79.57) > 0.1 ||
    validation$n_total != 154L
  ) {
    stop("Study 55 Welch reconstruction does not match the registered result.", call. = FALSE)
  }
  
  validation
}

study_55_prior_grid <- function(priors) {
  grid <- priors |>
    dplyr::filter(
      .data$prior_family == "welch_anova_normal",
      .data$param == "contrast_prior_sd"
    ) |>
    dplyr::transmute(
      prior_label = .data$prior_label,
      rscale = as.numeric(.data$value)
    ) |>
    dplyr::mutate(
      prior_order = match(.data$prior_label, c("narrow", "primary", "wide"))
    ) |>
    dplyr::arrange(.data$prior_order) |>
    dplyr::select(-"prior_order")
  
  if (
    nrow(grid) != 3L ||
    !setequal(grid$prior_label, c("narrow", "primary", "wide")) ||
    any(!is.finite(grid$rscale)) ||
    any(grid$rscale <= 0)
  ) {
    stop("Study 55 requires three valid welch_anova_normal priors.", call. = FALSE)
  }
  grid
}

study_55_formulas <- function() {
  list(model_null = brms::bf(outcome_z ~ 0 + Intercept,
                             sigma ~ 0 + condition),
       model_alt = brms::bf(outcome_z ~ 0 + Intercept + contrast_1 + contrast_2 + contrast_3,
                            sigma ~ 0 + condition)
       )
}

study_55_null_priors <- function() {
  c(
    brms::set_prior("normal(0, 1)", class = "b", coef = "Intercept"),
    brms::set_prior("normal(0, 1)", class = "b", dpar = "sigma")
  )
}

study_55_alt_priors <- function(rscale) {
  c(
    brms::set_prior("normal(0, 1)", class = "b", coef = "Intercept"),
    brms::set_prior(
      paste0("normal(0, ", format(rscale, digits = 17), ")"),
      class = "b",
      coef = c("contrast_1", "contrast_2", "contrast_3")
    ),
    brms::set_prior("normal(0, 1)", class = "b", dpar = "sigma")
  )
}

fit_study_55_null <- function(data, seed = 123) {
  formulas <- study_55_formulas()
  
  brms::brm(
    formula = formulas$model_null,
    data = data,
    family = stats::gaussian(),,
    prior = study_55_null_priors(),
    chains = 4,
    iter = 12000,
    warmup = 2000,
    cores = 4,
    seed = seed,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    save_pars = brms::save_pars(all = TRUE),
    refresh = 500
  )
}

fit_study_55_alt <- function(data, rscale, seed) {
  formulas <- study_55_formulas()
  
  brms::brm(
    formula = formulas$model_alt,
    data = data,
    family = stats::gaussian(),,
    prior = study_55_alt_priors(rscale),
    chains = 4,
    iter = 12000,
    warmup = 2000,
    cores = 4,
    seed = seed,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    save_pars = brms::save_pars(all = TRUE),
    refresh = 500
  )
}

validate_study_55_fit <- function(fit, model_label) {
  fit_summary <- summary(fit)
  fixed <- fit_summary$fixed
  
  if (any(!is.finite(fixed[, "Rhat"])) || max(fixed[, "Rhat"]) > 1.01) {
    stop("Study 55 ", model_label, " has an unacceptable R-hat.", call. = FALSE)
  }
  
  sampler <- brms::nuts_params(fit)
  divergences <- sum(sampler$Parameter == "divergent__" & sampler$Value == 1)
  
  if (divergences > 0L) {
    stop("Study 55 ", model_label, " has ", divergences, " divergences.", call. = FALSE)
  }
  
  invisible(TRUE)
}

bridge_study_55_pair <- function(
    null_fit,
    alt_fit,
    repetitions = 5,
    cores = 4) {
  
  estimates <- purrr::map_dfr(seq_len(repetitions), function(i) {
    null_bridge <- bridgesampling::bridge_sampler(null_fit, silent = TRUE, cores = cores)
    alt_bridge <- bridgesampling::bridge_sampler(alt_fit, silent = TRUE, cores = cores)
    
    log_bf10 <- as.numeric(alt_bridge$logml - null_bridge$logml)
    
    tibble::tibble(
      repetition = i,
      log_bf10 = log_bf10,
      log10_bf10 = log_bf10 / log(10)
    )
  })
  
  mean_log_bf10 <- mean(estimates$log_bf10)
  mean_log10_bf10 <- mean(estimates$log10_bf10)
  
  tibble::tibble(
    bf10 = exp(mean_log_bf10),
    log_bf10 = mean_log_bf10,
    log10_bf10 = mean_log10_bf10,
    bf_error = stats::sd(estimates$log10_bf10),
    bridge_span_log10 = diff(range(estimates$log10_bf10)),
    bridge_repetitions = repetitions
  )
}

run_study_55_bayes_factors <- function(
    priors_path = "config/priors_wave2.csv",
    cache_dir = "outputs/intermediate/study_55_models") {
  
  data <- load_study_55_data()
  validation <- validate_study_55_reconstruction(data)
  priors <- readr::read_csv(priors_path, show_col_types = FALSE)
  grid <- study_55_prior_grid(priors)
  
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  null_path <- file.path(cache_dir, "study_55_null.rds")
  
  if (file.exists(null_path)) {
    message("Study 55: reading cached null model.")
    null_fit <- readRDS(null_path)
  } else {
    message("Study 55: fitting null model.")
    null_fit <- fit_study_55_null(data)
    validate_study_55_fit(null_fit, "null model")
    saveRDS(null_fit, null_path)
  }
  
  purrr::map_dfr(seq_len(nrow(grid)), function(i) {
    prior_label <- grid$prior_label[[i]]
    rscale <- grid$rscale[[i]]
    alt_path <- file.path(cache_dir, paste0("study_55_alt_", prior_label, ".rds"))
    
    if (file.exists(alt_path)) {
      message("Study 55: reading cached ", prior_label, " alternative model.")
      alt_fit <- readRDS(alt_path)
    } else {
      message("Study 55: fitting ", prior_label, " alternative model.")
      alt_fit <- fit_study_55_alt(data, rscale = rscale, seed = 123 + i)
      validate_study_55_fit(alt_fit, paste0(prior_label, " alternative model"))
      saveRDS(alt_fit, alt_path)
    }
    
    validate_study_55_fit(null_fit, "null model")
    validate_study_55_fit(alt_fit, paste0(prior_label, " alternative model"))
    
    bridge <- bridge_study_55_pair(null_fit, alt_fit)
    
    tibble::tibble(
      claim_id = "study_55_claim_01",
      study_id = "study_55",
      prior_label = prior_label,
      rscale = rscale,
      bf10 = bridge$bf10,
      log_bf10 = bridge$log_bf10,
      log10_bf10 = bridge$log10_bf10,
      bf_error = bridge$bf_error,
      bridge_span_log10 = bridge$bridge_span_log10,
      bridge_repetitions = bridge$bridge_repetitions,
      model_null = "outcome_z ~ 0 + Intercept; sigma ~ 0 + condition",
      model_alt = "outcome_z ~ 0 + Intercept + contrast_1 + contrast_2 + contrast_3; sigma ~ 0 + condition",
      method = "heteroscedastic_gaussian_anova_bridge_sampling",
      f_reconstructed = validation$f_value,
      df1_reconstructed = validation$df1,
      df2_reconstructed = validation$df2,
      n_total = validation$n_total
    )
  })
}

compute_study_55_bayes_factors <- function(
    claim,
    priors,
    results_path = "outputs/intermediate/study_55_bayes_factors.csv") {
  
  if (!file.exists(results_path)) {
    stop("Study 55 cached results not found: ", results_path, ".", call. = FALSE)
  }
  
  cached <- readr::read_csv(results_path, show_col_types = FALSE) |>
    dplyr::filter(.data$claim_id == claim$claim_id) |>
    dplyr::mutate(
      prior_order = match(.data$prior_label, c("narrow", "primary", "wide"))
    ) |>
    dplyr::arrange(.data$prior_order)
  
  if (nrow(cached) != 3L || anyNA(cached$prior_order)) {
    stop("Study 55 requires three valid cached prior rows.", call. = FALSE)
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
      bf_error = row$bf_error,
      model_null = row$model_null,
      model_alt = row$model_alt,
      bf_family = "heteroscedastic_gaussian_anova",
      prior_family = "welch_anova_normal",
      method = row$method
    )
  })
}