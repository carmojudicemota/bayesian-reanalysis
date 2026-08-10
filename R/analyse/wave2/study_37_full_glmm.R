study_37_full_formula <- function() {
  brms::bf(
    Attrition ~
      h1_demo_age_nontrad1 +
      h1_demo_gender_male1 +
      h1_demo_ethn_Lat1 +
      h1_demo_ethn_Afr1 +
      h1_demo_parents_y_1 +
      h1_demo_addl_lang_EngOther1 +
      h1_demo_first_semest1 +
      h1_regents_avg_z +
      h1_ESSES_avg_z +
      h1_textbook_y_1 +
      h1_blackboard_y_1 +
      synchronous1 +
      HUMAN_attn +
      HUMAN_work +
      HUMAN_course +
      HUMAN_mental +
      HUMAN_sleep +
      HUMAN_family +
      (1 | section)
  )
}


study_37_null_formula <- function(claim_id) {
  switch(
    claim_id,
    study_37_claim_01 =
      brms::bf(
        Attrition ~
          h1_demo_age_nontrad1 +
          h1_demo_gender_male1 +
          h1_demo_ethn_Lat1 +
          h1_demo_ethn_Afr1 +
          h1_demo_parents_y_1 +
          h1_demo_addl_lang_EngOther1 +
          h1_demo_first_semest1 +
          h1_regents_avg_z +
          h1_ESSES_avg_z +
          h1_textbook_y_1 +
          h1_blackboard_y_1 +
          synchronous1 +
          HUMAN_attn +
          HUMAN_work +
          HUMAN_course +
          HUMAN_mental +
          HUMAN_sleep +
          (1 | section)
      ),
    
    study_37_claim_02 =
      brms::bf(
        Attrition ~
          h1_demo_age_nontrad1 +
          h1_demo_gender_male1 +
          h1_demo_ethn_Lat1 +
          h1_demo_ethn_Afr1 +
          h1_demo_parents_y_1 +
          h1_demo_addl_lang_EngOther1 +
          h1_demo_first_semest1 +
          h1_regents_avg_z +
          h1_ESSES_avg_z +
          h1_textbook_y_1 +
          h1_blackboard_y_1 +
          HUMAN_attn +
          HUMAN_work +
          HUMAN_course +
          HUMAN_mental +
          HUMAN_sleep +
          HUMAN_family +
          (1 | section)
      ),
    
    stop(
      "Unknown Study 37 claim: ",
      claim_id,
      call. = FALSE
    )
  )
}


load_study_37_full_glmm_data <- function(path = "data/raw/study_37/osf_data.csv") {
  if (!file.exists(path)) {
    stop(
      "Study 37 data file does not exist: ",
      path,
      call. = FALSE
    )
  }
  
  raw <- utils::read.csv(path, stringsAsFactors = FALSE)
  
  binary_predictors <- c(
    "h1_demo_age_nontrad1",
    "h1_demo_gender_male1",
    "h1_demo_ethn_Lat1",
    "h1_demo_ethn_Afr1",
    "h1_demo_parents_y_1",
    "h1_demo_addl_lang_EngOther1",
    "h1_demo_first_semest1",
    "h1_textbook_y_1",
    "h1_blackboard_y_1",
    "synchronous1",
    "HUMAN_attn",
    "HUMAN_work",
    "HUMAN_course",
    "HUMAN_mental",
    "HUMAN_sleep",
    "HUMAN_family"
  )
  
  required_columns <- c("Attrition",
                        "section",
                        "h1_regents_avg",
                        "h1_ESSES_avg",
                        binary_predictors
                        )
  
  missing_columns <- setdiff(required_columns,names(raw))
  
  if (length(missing_columns) > 0L) {
    stop(
      "Study 37 is missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  data <- raw[stats::complete.cases(
      raw[required_columns]
    ),
    required_columns,
    drop = FALSE
  ]
  
  data$Attrition <- as.integer(as.character(data$Attrition))
  
  if (!all(
    data$Attrition %in% c(0L, 1L)
  )) {
    stop(
      "Study 37 Attrition must contain only 0 and 1.",
      call. = FALSE
    )
  }
  
  for (variable in binary_predictors) {
    data[[variable]] <- factor(data[[variable]])
    
    if (nlevels(data[[variable]]) != 2L) {
      stop(
        "Study 37 predictor ",
        variable,
        " must have exactly two levels.",
        call. = FALSE
      )
    }
  }
  
  data$section <- factor(data$section)
  data$h1_regents_avg_z <- as.numeric(scale(as.numeric(data$h1_regents_avg)))
  data$h1_ESSES_avg_z <- as.numeric(scale(as.numeric(data$h1_ESSES_avg)))
  data$h1_regents_avg <- NULL
  data$h1_ESSES_avg <- NULL
  
  if (nlevels(data$section) != 10L) {
    warning(
      "Study 37 contains ",
      nlevels(data$section),
      " sections rather than the expected 10."
    )
  }
  data
}


study_37_focal_coefficient <- function(claim_id) {
  switch(
    claim_id,
    study_37_claim_01 = "HUMAN_family1",
    study_37_claim_02 = "synchronous11",
    stop("Unknown Study 37 claim: ", claim_id, call. = FALSE)
  )
}


study_37_model_priors <- function(claim_id,prior_values,include_focal,focal_prior_override = NULL) {

  focal_coefficient <- study_37_focal_coefficient(claim_id)
  nuisance_prior <- sprintf("normal(0, %.15g)",as.numeric(prior_values[["nuisance_sd"]]))
  intercept_prior <- sprintf("normal(0, %.15g)",as.numeric(prior_values[["intercept_sd"]]))
  random_prior <- sprintf("normal(0, %.15g)", as.numeric(prior_values[["random_sd"]]))
  focal_prior <- if (is.null(focal_prior_override)) {
    sprintf("normal(0, %.15g)",as.numeric(prior_values[["focal_sd"]]))
  } else {
    focal_prior_override
  }
  model_priors <- c(
    brms::set_prior(prior = nuisance_prior,class = "b"),
    brms::set_prior(prior = intercept_prior,class = "Intercept"),
    brms::set_prior(prior = random_prior,class = "sd")
    )
  if (isTRUE(include_focal)) {
    model_priors <- c(
      model_priors,
      brms::set_prior(prior = focal_prior,class = "b",coef = focal_coefficient)
      )
  }
  
  model_priors
}

bridge_study_37_pair <- function(model_pair,repetitions = 5L,cores = 1L) {
  purrr::map_dfr(seq_len(repetitions),function(repetition) {
      full_bridge <- brms::bridge_sampler(model_pair$full,silent = TRUE,cores = cores)
      null_bridge <- brms::bridge_sampler(model_pair$null,silent = TRUE,cores = cores)
      logml_full <- as.numeric(full_bridge$logml)
      logml_null <- as.numeric(null_bridge$logml)
      log_bf10 <- logml_full -logml_null
      
      tibble::tibble(
        repetition = repetition,
        logml_full = logml_full,
        logml_null = logml_null,
        log_bf10 = log_bf10,
        log10_bf10 = log_bf10 / log(10),
        bf10 = if (log_bf10 < log(.Machine$double.xmax)) {exp(log_bf10)} else {Inf}
        )
    }
  )
}


load_study_37_bayes_factor_results <- function(
    path = "outputs/tables/study_37_full_glmm_bayes_factors.csv") {
  if (!file.exists(path)) {
    stop("Study 37 Bayes-factor results do not exist: ",path,call. = FALSE)
  }
  
  results <- readr::read_csv(path,show_col_types = FALSE)
  required_columns <- c(
    "claim_id",
    "study_id",
    "prior_label",
    "rscale",
    "bf10",
    "log_bf10",
    "log10_bf10",
    "bf_error",
    "model_null",
    "model_alt",
    "method",
    "bridge_sd_log10",
    "bridge_span_log10",
    "bridge_repetitions"
  )
  missing_columns <- setdiff(required_columns,names(results))
  if (length(missing_columns) > 0L) {
    stop(
      "Study 37 results are missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  expected_claims <- c("study_37_claim_01","study_37_claim_02")
  expected_priors <- c("narrow","primary","wide")
  results <- results |>
    dplyr::filter(.data$claim_id %in% expected_claims)
  
  if (nrow(results) != 6L) {
    stop(
      "Study 37 should contain exactly six result rows.",
      call. = FALSE
    )
  }
  
  if (!setequal(unique(results$prior_label),expected_priors)) {
    stop(
      "Study 37 must contain narrow, primary and wide priors.",
      call. = FALSE
    )
  }
  
  if (any(!is.finite(results$log10_bf10))) {
    stop(
      "Study 37 contains non-finite log10 Bayes factors.",
      call. = FALSE
    )
  }
  
  if (any(results$bridge_span_log10 >= 0.10)) {
    stop(
      "At least one Study 37 bridge estimate is unstable.",
      call. = FALSE
    )
  }
  
  results
}


compute_study_37_bayes_factors <- function(claim,priors = NULL) {
  
  results <- load_study_37_bayes_factor_results()
  claim_results <- results |>
    dplyr::filter(.data$claim_id == claim$claim_id)
  
  if (nrow(claim_results) != 3L) {
    stop("Expected three Study 37 prior rows for ",claim$claim_id,".",call. = FALSE)
  }
  
  purrr::map_dfr(seq_len(nrow(claim_results)),
                 function(i) {
                   row <- claim_results[i,,drop = FALSE]
                   wave2_row(claim = claim,
                             prior_label = row$prior_label[[1]],
                             rscale = row$rscale[[1]],
                             bf10 = row$bf10[[1]],
                             log_bf10 = row$log_bf10[[1]],
                             log10_bf10 = row$log10_bf10[[1]],
                             bf_error = row$bf_error[[1]],
                             model_null = row$model_null[[1]],
                             model_alt = row$model_alt[[1]],
                             bf_family = "full_bayesian_glmm",
                             prior_family = "glmm_full",
                             method = "full_bayesian_glmm_bridge_sampling"
                             )
    }
  )
}

