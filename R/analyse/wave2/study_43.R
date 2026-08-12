study_43_contrast_name <- function(claim_id) {
  switch(claim_id, study_43_claim_01 = "testing_minus_restudy_reviewed",
    study_43_claim_02 = "reviewed_minus_unreviewed",
    stop("Unknown Study 43 claim: ", claim_id, call. = FALSE)
  )
}

load_study_43_wave2_data <- function(path = "data/raw/study_43/Datafile.sav") {
  if (!file.exists(path)) stop("Study 43 data file does not exist.", call. = FALSE)
  
  required <- c("Crit_Score_Testing_old", "Crit_Score_Testing_New",
                "Crit_Score_Restudy_old", "Crit_Score_Restudy_New",
                "Lecture", "questiontype_crit"
                )
  
  raw <- haven::read_sav(path)
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0L) stop("Study 43 is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  data <- raw |>
    dplyr::select(dplyr::all_of(required)) |>
    dplyr::mutate(
      dplyr::across(dplyr::starts_with("Crit_Score"), as.numeric),
      Lecture = factor(Lecture),
      questiontype_crit = factor(questiontype_crit)
    ) |>
    tidyr::drop_na() |>
    dplyr::mutate(
      testing_minus_restudy_reviewed = Crit_Score_Testing_old - Crit_Score_Restudy_old,
      reviewed_minus_unreviewed = (Crit_Score_Testing_old + Crit_Score_Restudy_old) / 2 - (Crit_Score_Testing_New + Crit_Score_Restudy_New) / 2
    )
  
  old_contrasts <- options("contrasts")
  on.exit(options(old_contrasts), add = TRUE)
  options(contrasts = c("contr.sum", "contr.poly"))
  nuisance <- stats::model.matrix(~ Lecture + questiontype_crit, data = data)[, -1, drop = FALSE]
  nuisance <- as.data.frame(nuisance)
  names(nuisance) <- paste0("nuisance_", seq_along(nuisance))
  
  dplyr::bind_cols(data, nuisance)
}

validate_study_43_claim <- function(data, claim_id) {
  outcome <- study_43_contrast_name(claim_id)
  formula <- stats::as.formula(paste0(outcome, " ~ Lecture + questiontype_crit"))
  old_contrasts <- options("contrasts")
  on.exit(options(old_contrasts), add = TRUE)
  options(contrasts = c("contr.sum", "contr.poly"))
  fit <- stats::lm(formula, data = data)
  coefficient <- summary(fit)$coefficients["(Intercept)", ]
  
  tibble::tibble(claim_id = claim_id,
                 estimate = unname(coefficient["Estimate"]),
                 se = unname(coefficient["Std. Error"]),
                 t_value = unname(coefficient["t value"]),
                 df = stats::df.residual(fit),
                 f_value = unname(coefficient["t value"])^2,
                 p_value = 2 * stats::pt(abs(unname(coefficient["t value"])), stats::df.residual(fit), lower.tail = FALSE)
                 )
}


study_43_model_data <- function(data, claim_id) {
  outcome_name <- study_43_contrast_name(claim_id)
  outcome <- data[[outcome_name]]
  outcome_sd <- stats::sd(outcome)
  
  if (!is.finite(outcome_sd) || outcome_sd <= 0) {
    stop("Study 43 contrast has invalid standard deviation.", call. = FALSE)
  }
  
  old_contrasts <- options("contrasts")
  on.exit(options(old_contrasts), add = TRUE)
  options(contrasts = c("contr.sum", "contr.poly"))
  
  design <- stats::model.matrix(~ Lecture + questiontype_crit,data = data)
  nuisance <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  nuisance <- as.data.frame(nuisance)
  names(nuisance) <- paste0("nuisance_", seq_len(ncol(nuisance)))
  
  model_data <- dplyr::bind_cols(
    tibble::tibble(outcome_z = outcome / outcome_sd,focal = 1),nuisance)
  list(data = model_data,nuisance_names = names(nuisance),outcome_sd = outcome_sd,original_outcome = outcome)
}

study_43_prior_grid <- function(priors) {
  grid <- priors |>
    dplyr::filter(.data$prior_family == "glmm_normal", .data$param == "prior_sd") |>
    dplyr::mutate(prior_order = match(.data$prior_label, c("narrow", "primary", "wide"))) |>
    dplyr::arrange(.data$prior_order)
  
  if (nrow(grid) != 3L || anyNA(grid$prior_order)) {
    stop("Study 43 requires narrow, primary and wide normal prior rows.", call. = FALSE)
  }
  
  grid
}

study_43_priors <- function(focal_sd, include_focal) {
  priors <- c(brms::set_prior("normal(0, 1)", class = "b"),
              brms::set_prior("student_t(3, 0, 1)", class = "sigma")
  )
  
  if (include_focal) {
    priors <- c(priors, brms::set_prior(sprintf("normal(0, %.15g)", focal_sd), class = "b", coef = "focal"))
  }
  
  priors
}

fit_study_43_models <- function(claim_id, prior_label, focal_sd, data = load_study_43_wave2_data(),
                                iter = 6000, warmup = 2000, chains = 4,
                                cores = min(chains, getOption("bayesian_reanalysis.cores", getOption("mc.cores", 1L))), seed = 123) {
  prepared <- study_43_model_data(data, claim_id)
  formulas <- study_43_formulas(prepared$nuisance_names)
  
  full <- brms::brm(formulas$full, data = prepared$data, family = gaussian(),
                    prior = study_43_priors(focal_sd, TRUE),
                    iter = iter, warmup = warmup, chains = chains, cores = cores, seed = seed,
                    save_pars = brms::save_pars(all = TRUE), refresh = 0
  )
  
  null <- brms::brm(formulas$null, data = prepared$data, family = gaussian(),
                    prior = study_43_priors(focal_sd, FALSE),
                    iter = iter, warmup = warmup, chains = chains, cores = cores, seed = seed + 1,
                    save_pars = brms::save_pars(all = TRUE), refresh = 0
  )
  
  list(full = full, null = null, model_full = deparse(formulas$full$formula), model_null = deparse(formulas$null$formula))
}

bridge_study_43_models <- function(models, repetitions = 5L,
                                   cores = getOption("bayesian_reanalysis.cores", getOption("mc.cores", 1L))) {
  estimates <- purrr::map_dfr(seq_len(repetitions), function(i) {
    full_bridge <- bridgesampling::bridge_sampler(models$full, silent = TRUE, cores = cores)
    null_bridge <- bridgesampling::bridge_sampler(models$null, silent = TRUE, cores = cores)
    log_bf10 <- as.numeric(full_bridge$logml - null_bridge$logml)
    
    tibble::tibble(repetition = i,
                   log_bf10 = log_bf10,
                   log10_bf10 = log_bf10 / log(10),
                   bf10 = exp(log_bf10)
                   )
  })
  
  tibble::tibble(bf10 = exp(mean(estimates$log_bf10)),
                 log_bf10 = mean(estimates$log_bf10),
                 log10_bf10 = mean(estimates$log10_bf10),
                 bf_error = stats::sd(estimates$log10_bf10),
                 bridge_span_log10 = diff(range(estimates$log10_bf10)),
                 bridge_repetitions = repetitions
                 )
}

run_study_43_bayes_factors <- function(
    priors_path = "config/priors_wave2.csv",
    output_path = "outputs/tables/study_43_bayes_factors.csv") {
  
  priors <- readr::read_csv(priors_path, show_col_types = FALSE)
  grid <- study_43_prior_grid(priors)
  data <- load_study_43_wave2_data()
  results <- purrr::map_dfr(c("study_43_claim_01", "study_43_claim_02"), function(claim_id) {
    purrr::map_dfr(seq_len(nrow(grid)), function(i) {
      models <- fit_study_43_models(claim_id = claim_id,
                                    prior_label = grid$prior_label[[i]],
                                    focal_sd = grid$value[[i]],
                                    data = data)
      result <- bridge_study_43_models(models)
      dplyr::mutate(result,
                    claim_id = claim_id,
                    study_id = "study_43",
                    prior_label = grid$prior_label[[i]],
                    rscale = grid$value[[i]],
                    model_null = paste(models$model_null, collapse = ""),
                    model_alt = paste(models$model_full, collapse = ""),
                    method = "adjusted_contrast_gaussian_bridge_sampling",
                    .before = 1
                    )
    })
  })
  
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(results, output_path)
  invisible(results)
}

compute_study_43_bayes_factors <- function(claim,priors,results_path = "outputs/intermediate/study_43_bayes_factors.csv") {
  
  if (!file.exists(results_path)) {
    stop("Study 43 cached results not found: ",results_path,call. = FALSE)
  }
  
  cached <- readr::read_csv(results_path,show_col_types = FALSE
  ) |>
    dplyr::filter(.data$claim_id == claim$claim_id) |>
    dplyr::mutate(prior_order = match(.data$prior_label,c("narrow", "primary", "wide"))) |>
    dplyr::arrange(.data$prior_order)
  
  if (nrow(cached) != 3L ||anyNA(cached$prior_order)) {
    stop("Study 43 requires three valid cached prior rows for ",claim$claim_id,".",call. = FALSE)
  }
  
  purrr::pmap_dfr(cached,function(...) {
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
        bf_family = "gaussian_adjusted_contrast",
        prior_family = "normal_focal_bridge",
        method = row$method
      )
    }
  )
}

