load_study_10_wave2_data <- function(
    outcome_column,
    path = "data/raw/study_10/master.anonymizedOSF.csv") {
  
  if (!file.exists(path)) {
    stop("Study 10 data file does not exist: ", path, call. = FALSE)
  }
  
  allowed_outcomes <- c("GHQc.total", "QHQb.total")
  
  if (!outcome_column %in% allowed_outcomes) {
    stop(
      "Unsupported Study 10 outcome: ", outcome_column,
      ". Expected one of: ", paste(allowed_outcomes, collapse = ", "),
      call. = FALSE
    )
  }
  
  raw <- readr::read_csv(
    path,
    show_col_types = FALSE,
    name_repair = "unique_quiet"
  )
  
  required_columns <- c("university", "replicate", outcome_column)
  missing_columns <- setdiff(required_columns, names(raw))
  
  if (length(missing_columns) > 0L) {
    stop(
      "Study 10 is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  data <- raw |>
    dplyr::filter(.data$university %in% c(1, 3)) |>
    dplyr::transmute(
      outcome = as.numeric(.data[[outcome_column]]),
      A = factor(.data$replicate, levels = c("a", "b", "c")),
      B = factor(
        .data$university,
        levels = c(1, 3),
        labels = c("university_1", "university_3")
      )
    ) |>
    tidyr::drop_na() |>
    as.data.frame()
  
  if (nrow(data) != 1049L) {
    stop(
      "Unexpected Study 10 sample size for ",
      outcome_column, ": ", nrow(data), ".",
      call. = FALSE
    )
  }
  
  data
}

study_10_outcome <- function(claim_id) {
  switch(
    claim_id,
    study_10_claim_01 = "GHQc.total",
    study_10_claim_02 = "QHQb.total",
    stop("Unknown Study 10 claim: ", claim_id, call. = FALSE)
  )
}

study_10_expected_f <- function(claim_id) {
  switch(
    claim_id,
    study_10_claim_01 = 103.499,
    study_10_claim_02 = 25.505,
    stop("Unknown Study 10 claim: ", claim_id, call. = FALSE)
  )
}

extract_anova_value <- function(table, term, statistic) {
  row_names <- trimws(rownames(table))
  column_names <- trimws(colnames(table))
  
  row_index <- match(term, row_names)
  column_index <- match(statistic, column_names)
  
  if (is.na(row_index) || is.na(column_index)) {
    stop(
      "Could not extract term '", term,
      "' and statistic '", statistic,
      "'. Available rows: ", paste(row_names, collapse = " | "),
      ". Available columns: ", paste(column_names, collapse = " | "),
      call. = FALSE
    )
  }
  
  value <- as.numeric(table[row_index, column_index])
  
  if (length(value) != 1L || !is.finite(value)) {
    stop(
      "Extracted ANOVA value is not finite for term '",
      term, "' and statistic '", statistic, "'.",
      call. = FALSE
    )
  }
  
  value
}

validate_study_10_original_anova <- function(claim_id, data) {
  type_1_fit <- stats::aov(outcome ~ A * B, data = data)
  type_1_table <- summary(type_1_fit)[[1]]
  
  contrast_fit <- stats::lm(
    outcome ~ A * B,
    data = data,
    contrasts = list(A = "contr.sum", B = "contr.sum")
  )
  
  type_2_table <- car::Anova(contrast_fit, type = 2)
  type_3_table <- car::Anova(contrast_fit, type = 3)
  
  intercept_fit <- stats::lm(outcome ~ 1, data = data)
  semester_fit <- stats::lm(outcome ~ A, data = data)
  intercept_comparison <- stats::anova(intercept_fit, semester_fit)
  
  f_type_1 <- extract_anova_value(type_1_table, "A", "F value")
  f_type_2 <- extract_anova_value(type_2_table, "A", "F value")
  f_type_3 <- extract_anova_value(type_3_table, "A", "F value")
  ss_type_1 <- extract_anova_value(type_1_table, "A", "Sum Sq")
  
  f_intercept_vs_a <- as.numeric(intercept_comparison[2, "F"])
  ss_intercept_vs_a <- as.numeric(intercept_comparison[2, "Sum of Sq"])
  
  values_to_check <- c(
    f_type_1,
    f_type_2,
    f_type_3,
    f_intercept_vs_a,
    ss_type_1,
    ss_intercept_vs_a
  )
  
  if (any(!is.finite(values_to_check))) {
    stop(
      "Study 10 validation produced one or more non-finite values.",
      call. = FALSE
    )
  }
  
  result <- tibble::tibble(
    claim_id = claim_id,
    outcome = study_10_outcome(claim_id),
    n_total = nrow(data),
    residual_df = stats::df.residual(type_1_fit),
    f_type_1 = f_type_1,
    f_type_2 = f_type_2,
    f_type_3 = f_type_3,
    f_intercept_vs_a = f_intercept_vs_a,
    ss_type_1 = ss_type_1,
    ss_intercept_vs_a = ss_intercept_vs_a
  )
  
  expected_f <- study_10_expected_f(claim_id)
  
  if (abs(result$f_type_1 - expected_f) > 0.02) {
    stop(
      "Study 10 Type I validation failed for ",
      claim_id, ": expected approximately ", expected_f,
      ", obtained ", result$f_type_1, ".",
      call. = FALSE
    )
  }
  
  if (abs(result$ss_type_1 - result$ss_intercept_vs_a) > 1e-8) {
    stop(
      "Study 10 sequential numerator validation failed for ",
      claim_id, ".",
      call. = FALSE
    )
  }
  
  result
}

compute_study_10_bayes_factors <- function(claim, priors) {
  outcome_column <- study_10_outcome(claim$claim_id)
  data <- load_study_10_wave2_data(outcome_column)
  
  validate_study_10_original_anova(
    claim_id = claim$claim_id,
    data = data
  )
  
  compute_factorial_additive_main_effect_bfs(
    claim = claim,
    data = data,
    priors = priors
  )
}

compute_study_10_marginality_sensitivity <- function(claim, priors) {
  outcome_column <- study_10_outcome(claim$claim_id)
  data <- load_study_10_wave2_data(outcome_column)
  
  validate_study_10_original_anova(
    claim_id = claim$claim_id,
    data = data
  )
  
  compute_factorial_main_effect_bfs(
    claim = claim,
    data = data,
    priors = priors
  ) |>
    dplyr::mutate(
      analysis_role = "type3_coefficient_block_sensitivity"
    )
}