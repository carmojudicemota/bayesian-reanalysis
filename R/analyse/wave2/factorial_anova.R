find_factorial_model <- function(bayes_factors,model_name,claim_id) {
  model_table <- BayesFactor::extractBF(bayes_factors)
  available_models <- rownames(model_table)
  
  if (!model_name %in% available_models) {
    stop(
      "Model '",model_name,"' was not returned for ",claim_id,
      ". Available models: ",paste(available_models,collapse = " | "),
      call. = FALSE
    )
  }
  
  bayes_factors[model_name]
}

factorial_prior_grid <- function(priors) {
  prior_grid <- priors |>
    dplyr::filter(
      .data$prior_family == "anova_cauchy",
      .data$param == "rscale_fixed"
    ) |>
    dplyr::mutate(
      value = suppressWarnings(as.numeric(.data$value)),
      prior_order = match(
        .data$prior_label,
        c("narrow","primary","wide")
      )
    ) |>
    dplyr::arrange(.data$prior_order)
  
  if (
    nrow(prior_grid) != 3L ||
    anyNA(prior_grid$prior_order) ||
    any(!is.finite(prior_grid$value)) ||
    any(prior_grid$value <= 0)
  ) {
    stop(
      "The factorial analysis requires exactly the narrow, primary and wide ANOVA priors.",
      call. = FALSE
    )
  }
  
  prior_grid
}

compute_factorial_model_pair_bfs <- function(
    claim,
    data,
    priors,
    formula,
    model_null,
    model_alt) {
  
  data <- as.data.frame(data)
  required_columns <- c("outcome", "A", "B")
  missing_columns <- setdiff(required_columns, names(data))
  
  if (length(missing_columns) > 0L) {
    stop(
      claim$claim_id, " factorial data are missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  prior_grid <- factorial_prior_grid(priors)
  
  purrr::map_dfr(seq_len(nrow(prior_grid)), function(i) {
    set.seed(123 + i)
    prior_label <- prior_grid$prior_label[[i]]
    rscale <- as.numeric(prior_grid$value[[i]])
    all_models <- BayesFactor::generalTestBF(
      formula = formula,
      data = data,
      whichModels = "all",
      rscaleFixed = rscale,
      progress = FALSE
    )
    
    alternative_bf <- find_factorial_model(
      bayes_factors = all_models,
      model_name = model_alt,
      claim_id = claim$claim_id
    )
    
    null_bf <- find_factorial_model(
      bayes_factors = all_models,
      model_name = model_null,
      claim_id = claim$claim_id
    )
    
    comparison <- alternative_bf / null_bf
    extracted <- BayesFactor::extractBF(comparison)
    
    wave2_row(
      claim = claim,
      prior_label = prior_label,
      rscale = rscale,
      bf10 = extracted$bf[[1]],
      bf_error = extracted$error[[1]],
      prior_family = "anova_cauchy",
      model_null = paste0("outcome ~ ", model_null),
      model_alt = paste0("outcome ~ ", model_alt)
    )
  })
}

compute_factorial_main_effect_bfs <- function(claim, data, priors) {
  compute_factorial_model_pair_bfs(
    claim = claim,
    data = data,
    priors = priors,
    formula = outcome ~ A + B + A:B,
    model_null = "B + A:B",
    model_alt = "A + B + A:B"
  )
}

compute_factorial_additive_main_effect_bfs <- function(claim, data, priors) {
  compute_factorial_model_pair_bfs(
    claim = claim,
    data = data,
    priors = priors,
    formula = outcome ~ A + B,
    model_null = "B",
    model_alt = "A + B"
  )
}




