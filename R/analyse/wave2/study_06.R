load_study_06_wave2_data <- function(path = paste0("data/raw/study_06/","Social_Annotation_and_SOB_SOC_Data.sav")) {
  
  candidate_paths <- c(path, paste0("data/raw/study_06/","Social Annotation and SOB_SOC Data.sav"))
  
  existing_paths <- candidate_paths[file.exists(candidate_paths)]
  if (length(existing_paths) == 0L) {stop("Study 06 data file does not exist.",call. = FALSE)}
  
  raw <- haven::read_sav(existing_paths[[1]])
  required_columns <- c("BCBSTotalScore","Condition","Course")
  missing_columns <- setdiff(required_columns,names(raw))
  
  if (length(missing_columns) > 0L) {
    stop(
      "Study 06 is missing required columns: ",
      paste(missing_columns,collapse = ", "),
      call. = FALSE
    )
  }
  
  data <-
    tibble::tibble(outcome =as.numeric(raw$BCBSTotalScore),
                   condition_numeric = as.numeric(raw$Condition),
                   course_numeric = as.numeric(raw$Course),
                   condition_label = as.character(haven::as_factor(raw$Condition)),
                   course_label = as.character(haven::as_factor(raw$Course))
                   ) |>
    tidyr::drop_na(outcome,condition_numeric,course_numeric)
  
  if (dplyr::n_distinct(data$condition_numeric) != 2L) {
    stop("Study 06 Condition must have exactly two levels.",call. = FALSE)
  }
  if (dplyr::n_distinct(data$course_numeric) != 2L) {
    stop("Study 06 Course must have exactly two levels.",call. = FALSE)
  }
  
  condition_key <- data |> dplyr::distinct(condition_numeric,condition_label)
  individual_value <- condition_key |>
    dplyr::filter(
      grepl("individual|control",condition_label, ignore.case = TRUE)) |>
    dplyr::pull(condition_numeric)
  social_value <- condition_key |>
    dplyr::filter(
      grepl("social|perusall|annotation",condition_label,ignore.case = TRUE),
      !grepl("individual",condition_label,ignore.case = TRUE)
    ) |>
    dplyr::pull(condition_numeric)
  
  if (length(individual_value) != 1L || length(social_value) != 1L) {
    stop("Study 06 annotation conditions could not be identified uniquely.",call. = FALSE)
  }
  
  data <- data |>
    dplyr::mutate(
      condition_social = as.numeric(condition_numeric == social_value),
      course_other = as.numeric(course_numeric != 1),
      condition_social_course_other = condition_social * course_other,
      condition = factor(condition_social, 
                         levels = c(0,1),
                         labels = c("individual","social")),
      course = factor(course_other,
                      levels = c(0,1),
                      labels = c("psyc111","other_course"))
    ) |>
    dplyr::select(outcome,
                  condition,
                  course,
                  condition_social,
                  course_other,
                  condition_social_course_other
    ) |>
    as.data.frame()
  
  cell_counts <- table(data$condition,data$course)
  
  if (any(cell_counts == 0L)) {
    stop("Study 06 contains an empty factorial cell.",call. = FALSE)
  }
  
  data
}

compute_study_06_bayes_factors <- function(claim,priors) {
  
  data <- load_study_06_wave2_data()
  prior_grid <- priors |>
    dplyr::filter(.data$prior_family == "regression_mixture_g",
                  .data$param == "rscale_cont"
                  ) |>
    dplyr::mutate(prior_order = match(.data$prior_label,
          c("narrow", "primary", "wide"))
    ) |>
    dplyr::arrange(.data$prior_order)
  
  if (nrow(prior_grid) != 3L || anyNA(prior_grid$prior_order)
  ) {
    stop(
      "Study 06 requires narrow, primary and wide regression prior rows.",
      call. = FALSE
    )
  }
  
  model_null <-
    outcome ~
    course_other +
    condition_social_course_other
  
  model_alt <-
    outcome ~
    condition_social +
    course_other +
    condition_social_course_other
  
  purrr::map_dfr(
    seq_len(nrow(prior_grid)),
    function(i) {
      set.seed(123 + i)
      prior_label <-prior_grid$prior_label[[i]]
      fixed_scale <- as.numeric(prior_grid$value[[i]])
      null_bf <- BayesFactor::lmBF(formula = model_null,
                                   data = data,
                                   rscaleCont = fixed_scale,
                                   progress =FALSE
                                   )
      alt_bf <- BayesFactor::lmBF(formula = model_alt,
                                  data = data,
                                  rscaleCont = fixed_scale,
                                  progress = FALSE
                                  )
      comparison <- alt_bf / null_bf
      extracted <- BayesFactor::extractBF(comparison)
      wave2_row(claim =claim,
                prior_label = prior_label,
                rscale =fixed_scale,
                bf10 = extracted$bf[[1]],
                bf_error = extracted$error[[1]],
                model_null = paste(deparse(model_null),collapse = ""),
                model_alt = paste(deparse(model_alt),collapse = ""),
                bf_family = "gaussian_regression_simple_effect",
                prior_family = "regression_mixture_g",
                method = "full_factorial_gaussian_simple_effect"
      )
    }
  )
}

