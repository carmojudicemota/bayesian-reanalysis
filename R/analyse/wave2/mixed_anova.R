compute_mixed_anova_model_pair <- function(
    claim,
    data,
    model_null,
    model_alt,
    priors,
    participant_column = "participant") {
  
  data <- as.data.frame(data)
  
  required_columns <- unique(
    c(all.vars(model_null),all.vars(model_alt),participant_column)
  )
  
  missing_columns <- setdiff(required_columns,names(data))
  
  if (length(missing_columns) > 0L) {
    stop(
      claim$claim_id,
      " mixed-ANOVA data are missing: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  if (!is.factor(data[[participant_column]])) {
    data[[participant_column]] <- factor(
      data[[participant_column]]
    )
  }
  
  fixed_priors <- priors |>
    dplyr::filter(
      .data$prior_family == "anova_cauchy",
      .data$param == "rscale_fixed"
    ) |>
    dplyr::mutate(
      prior_order = match(
        .data$prior_label,
        c("narrow","primary","wide")
      )
    ) |>
    dplyr::arrange(.data$prior_order)
  
  if (
    nrow(fixed_priors) != 3L ||
    anyNA(fixed_priors$prior_order)
  ) {
    stop(
      "The mixed ANOVA requires exactly the narrow, ",
      "primary and wide fixed-effect priors.",
      call. = FALSE
    )
  }
  
  random_scale <- priors |>
    dplyr::filter(
      .data$prior_family == "anova_cauchy",
      .data$param == "rscale_random",
      .data$prior_label == "primary"
    ) |>
    dplyr::pull(
      .data$value
    )
  
  if (length(random_scale) != 1L) {
    stop(
      "The mixed ANOVA requires exactly one primary ",
      "rscale_random prior.",
      call. = FALSE
    )
  }
  
  purrr::map_dfr(
    seq_len(nrow(fixed_priors)),
    function(i) {
      set.seed(123 + i)
      prior_label <- fixed_priors$prior_label[i]
      fixed_scale <- as.numeric(fixed_priors$value[i])
      null_bf <- BayesFactor::lmBF(formula = model_null,
                                   data = data,
                                   whichRandom = participant_column,
                                   rscaleFixed = fixed_scale,
                                   rscaleRandom = random_scale,
                                   progress = FALSE
                                   )
      
      alternative_bf <- BayesFactor::lmBF(formula = model_alt,
                                          data = data,
                                          whichRandom = participant_column,
                                          rscaleFixed = fixed_scale,
                                          rscaleRandom = random_scale,
                                          progress = FALSE
                                          )
      
      comparison <- alternative_bf / null_bf
      extracted <- BayesFactor::extractBF(comparison)
      wave2_row(
        claim = claim,
        prior_label = prior_label,
        rscale = fixed_scale,
        bf10 = extracted$bf[1],
        bf_error = extracted$error[1],
        model_null = paste(
          deparse(model_null),
          collapse = ""
        ),
        model_alt = paste(
          deparse(model_alt),
          collapse = ""
        )
      )
    }
  )
}