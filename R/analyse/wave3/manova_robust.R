source("R/analyse/wave3/wave3_helpers.R")
for (f in c("R/analyse/wave3/study_03.R", "R/analyse/wave3/study_18.R", "R/analyse/wave3/study_44.R")) {
  if (file.exists(f)) source(f)
}


manova_robust_spec <- function(claim_id) {
  switch(
    claim_id,
    study_03_claim_01 = list(
      load = load_study_03_data, outcomes = study_03_outcomes(),
      full = "syllabus_snapshot * instructor_snapshot",
      null = "syllabus_snapshot + instructor_snapshot"
    ),
    study_03_claim_02 = list(
      load = load_study_03_data, outcomes = study_03_outcomes(),
      full = "syllabus_snapshot * instructor_snapshot",
      null = "instructor_snapshot + syllabus_snapshot:instructor_snapshot"
    ),
    study_18_claim_01 = list(
      load = load_study_18_data, outcomes = paste0("skill_", 1:8),
      full = "subject_group",
      null = "1"
    ),
    study_44_claim_01 = list(
      load = load_study_44_data, outcomes = study_44_outcomes(),
      full = "participant_status * rating_magnitude * quizzing_comment",
      null = "participant_status * rating_magnitude * quizzing_comment - rating_magnitude"
    ),
    study_44_claim_02 = list(
      load = load_study_44_data, outcomes = study_44_outcomes(),
      full = "participant_status * rating_magnitude * quizzing_comment",
      null = "participant_status * rating_magnitude * quizzing_comment - quizzing_comment"
    ),
    stop("No Student-t robustness spec for ", claim_id, call. = FALSE)
  )
}

manova_robust_claim_ids <- function() {
  c("study_03_claim_01", "study_03_claim_02", "study_18_claim_01",
    "study_44_claim_01", "study_44_claim_02")
}


manova_robust_standardize <- function(data, outcomes) {
  for (o in outcomes) data[[o]] <- as.numeric(scale(data[[o]]))
  data
}

fit_manova_robust <- function(data, outcomes, rhs, focal_scale = 0.5, seed = 123,
                              family = brms::student(),
                              chains = 4, iter = 4000, warmup = 1000, cores = 4) {
  form <- stats::as.formula(paste0("mvbind(", paste(outcomes, collapse = ", "), ") ~ ", rhs))
  model <- brms::bf(form) + brms::set_rescor(TRUE)
  priors <- brms::get_prior(model, data = data, family = family)
  priors$prior[priors$class == "Intercept" & priors$coef == ""] <- "student_t(3, 0, 2.5)"
  priors$prior[priors$class == "b" & priors$coef == ""] <- paste0("normal(0, ", focal_scale, ")")
  priors$prior[priors$class == "rescor"] <- "lkj(1)"
  brms::brm(
    formula = model,
    data = data,
    family = family,
    prior = priors,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    seed = seed,
    save_pars = brms::save_pars(all = TRUE),
    control = list(adapt_delta = 0.95),
    refresh = 0
  )
}

bridge_manova_robust <- function(full_fit, null_fit, repetitions = 5L, cores = 1L) {
  estimates <- vapply(seq_len(repetitions), function(i) {
    full_bridge <- brms::bridge_sampler(full_fit, silent = TRUE, cores = cores)
    null_bridge <- brms::bridge_sampler(null_fit, silent = TRUE, cores = cores)
    as.numeric(full_bridge$logml - null_bridge$logml)
  }, numeric(1))
  log_bf10 <- stats::median(estimates)
  list(
    log_bf10 = log_bf10,
    log10_bf10 = log_bf10 / log(10),
    bf10 = exp(log_bf10),
    bridge_sd_log10 = stats::sd(estimates / log(10)),
    bridge_span_log10 = diff(range(estimates / log(10))),
    bridge_repetitions = repetitions
  )
}

manova_robust_cache_path <- function() {"outputs/intermediate/manova_robust.csv"}
run_manova_robust <- function(claim_id, focal_scale = 0.5, repetitions = 5L, cache = manova_robust_cache_path()) {
  spec <- manova_robust_spec(claim_id)
  data <- manova_robust_standardize(spec$load(), spec$outcomes)
  full <- fit_manova_robust(data, spec$outcomes, spec$full, focal_scale, seed = 123)
  null <- fit_manova_robust(data, spec$outcomes, spec$null, focal_scale, seed = 124)
  br <- bridge_manova_robust(full, null, repetitions)
  row <- data.frame(
    claim_id = claim_id,
    study_id = sub("_claim.*", "", claim_id),
    focal_scale = focal_scale,
    bf10 = br$bf10,
    log10_bf10 = br$log10_bf10,
    bridge_sd_log10 = br$bridge_sd_log10,
    bridge_span_log10 = br$bridge_span_log10,
    bridge_repetitions = br$bridge_repetitions,
    method = "brms_manova_student_t_bridge",
    generated_at = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
  dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
  existing <- if (file.exists(cache)) utils::read.csv(cache, stringsAsFactors = FALSE) else NULL
  if (!is.null(existing)) existing <- existing[existing$claim_id != claim_id, , drop = FALSE]
  combined <- if (is.null(existing) || nrow(existing) == 0) row else rbind(existing, row)
  utils::write.csv(combined, cache, row.names = FALSE)
  invisible(row)
}


run_all_manova_robust <- function(claim_ids = manova_robust_claim_ids(),
                                  focal_scale = 0.5, repetitions = 5L) {
  purrr::map_dfr(claim_ids, function(cid) run_manova_robust(cid, focal_scale, repetitions))
}

diagnose_manova_robust <- function(claim_id, focal_scale = 0.5, repetitions = 5L) {
  spec <- manova_robust_spec(claim_id)
  data <- manova_robust_standardize(spec$load(), spec$outcomes)
  fit_pair <- function(family) {
    full <- fit_manova_robust(data, spec$outcomes, spec$full, focal_scale, seed = 123, family = family)
    null <- fit_manova_robust(data, spec$outcomes, spec$null, focal_scale, seed = 124, family = family)
    bridge_manova_robust(full, null, repetitions)$log10_bf10
  }
  tibble::tibble(
    claim_id = claim_id,
    focal_scale = focal_scale,
    gaussian_log10_bf10 = fit_pair(stats::gaussian()),
    student_t_log10_bf10 = fit_pair(brms::student())
  )
}

summarise_manova_robust <- function(primary = "outputs/tables/bayes_factor_results.csv",
                                    cache = manova_robust_cache_path(),
                                    out = "outputs/tables/manova_robust_sensitivity.csv",
                                    threshold = log10(3)) {
  if (!file.exists(cache)) {stop("No Student-t robustness cache at ", cache, "; run run_all_manova_robust() first.", call. = FALSE)}
  rob <- utils::read.csv(cache, stringsAsFactors = FALSE)
  cat_of <- function(x) ifelse(x > threshold, "H1", ifelse(x < -threshold, "H0", "Inconclusive"))
  pri <- readr::read_csv(primary, show_col_types = FALSE) |>
    dplyr::filter(.data$claim_id %in% rob$claim_id, .data$prior_label == "primary") |>
    dplyr::transmute(.data$claim_id, afbf_log10_bf10 = .data$log10_bf10)
  res <- pri |>
    dplyr::left_join(
      dplyr::transmute(rob, .data$claim_id,
                       student_t_log10_bf10 = .data$log10_bf10,
                       bridge_sd_log10 = .data$bridge_sd_log10),
      by = "claim_id"
    ) |>
    dplyr::mutate(
      afbf_category = cat_of(.data$afbf_log10_bf10),
      student_t_category = cat_of(.data$student_t_log10_bf10),
      category_agrees = .data$afbf_category == .data$student_t_category
    )
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(res, out)
  message("Student-t robustness: ", sum(res$category_agrees, na.rm = TRUE), "/",
          nrow(res), " claims agree with the AFBF category.")
  invisible(res)
}
