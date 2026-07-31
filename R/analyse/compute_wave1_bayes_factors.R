library(readr)
library(dplyr)
library(purrr)
library(BayesFactor)

compute_wave1_bayes_factors <- function(
    claims_path = "data/derived/claims.csv",
    priors_path = "config/priors_wave1.csv",
    output_path = "outputs/tables/bayes_factor_results_wave1.csv"
) {
  claims <- read_csv(claims_path, show_col_types = FALSE)
  priors <- read_csv(priors_path, show_col_types = FALSE)

  if (!"in_scope" %in% names(claims)) {
    claims$in_scope <- !is.na(claims$ali_result_status) &
      trimws(claims$ali_result_status) == "Fully reproducible"
  }

  allowed_sidedness <- c("one_sided", "two_sided", "omnibus", NA_character_)
  wave1_tests <- c(
    "one_sample_t_test", "paired_t_test",
    "independent_t_test", "repeated_measures_anova",
    "pearson_correlation"
  )

  wave1 <- claims |>
    filter(status == "ready", in_scope) |>
    filter(frequentist_test %in% wave1_tests) |>
    mutate(
      bf_family = if_else(frequentist_test == "pearson_correlation", "correlation", "t_test"),
      design = case_when(
        frequentist_test == "independent_t_test" ~ "independent",
        frequentist_test == "pearson_correlation" ~ "correlation",
        TRUE ~ "one_sample"
      ),
      t_for_bf = if_else(
        frequentist_test == "repeated_measures_anova" & is.na(t_value) & !is.na(f_value),
        sqrt(f_value), t_value
      ),
      df_for_bf = if_else(
        frequentist_test == "repeated_measures_anova" & is.na(t_df) & !is.na(f_df2),
        f_df2, t_df
      ),
      bf_sidedness = if_else(
        !is.na(p_sidedness) & p_sidedness == "one_sided",
        "one_sided", "two_sided"
      ),
      bf_direction = if_else(bf_sidedness == "one_sided", direction, NA_character_),
      observed_sign = case_when(
        bf_family == "correlation" & !is.na(r_value) & r_value < 0 ~ "negative",
        bf_family == "correlation" & !is.na(r_value) & r_value > 0 ~ "positive",
        bf_family == "correlation" & !is.na(r_value)              ~ "zero",
        !is.na(t_for_bf) & t_for_bf < 0                           ~ "negative",
        !is.na(t_for_bf) & t_for_bf > 0                           ~ "positive",
        !is.na(t_for_bf)                                          ~ "zero",
        TRUE                                                      ~ NA_character_
      ),
      direction_matches_observed = if_else(
        bf_sidedness == "one_sided" & !is.na(bf_direction) & !is.na(observed_sign),
        bf_direction == observed_sign,
        NA
      )
    )

  bad_label <- !(wave1$p_sidedness %in% allowed_sidedness)
  if (any(bad_label)) {
    stop(
      "Unrecognised p_sidedness label (expected one_sided/two_sided/omnibus/NA) for: ",
      paste(wave1$claim_id[bad_label], " = '", wave1$p_sidedness[bad_label], "'",
            sep = "", collapse = ", "),
      call. = FALSE
    )
  }

  undirected <- wave1$bf_sidedness == "one_sided" &
    !(wave1$bf_direction %in% c("positive", "negative"))
  if (any(undirected)) {
    stop(
      "One-sided claims with no recorded direction. Set `direction` to ",
      "'positive' or 'negative' in config/claim_map.csv for: ",
      paste(wave1$claim_id[undirected], collapse = ", "),
      call. = FALSE
    )
  }

  mismatched <- which(!is.na(wave1$direction_matches_observed) &
                        !wave1$direction_matches_observed)
  if (length(mismatched) > 0) {
    warning(
      "Recorded hypothesis direction disagrees with the recomputed sign for: ",
      paste(wave1$claim_id[mismatched], collapse = ", "),
      ". The recorded direction is used; record the discrepancy in the claim note.",
      call. = FALSE
    )
  }

  near_equal <- function(a, b) !is.na(a) & !is.na(b) & abs(a - b) < 1e-6

  rm_rows <- wave1$frequentist_test == "repeated_measures_anova"
  bad_rm <- rm_rows & !near_equal(wave1$f_df1, 1)
  if (any(bad_rm)) {
    stop(
      "repeated_measures_anova rows must have f_df1 == 1 to be reduced via sqrt(F). ",
      "Failing (missing or != 1): ",
      paste(wave1$claim_id[bad_rm], collapse = ", "),
      call. = FALSE
    )
  }

  ind_rows <- wave1$frequentist_test == "independent_t_test"
  bad_ind <- ind_rows & (is.na(wave1$n1) | is.na(wave1$n2) |
                           !near_equal(wave1$t_df, wave1$n1 + wave1$n2 - 2))
  if (any(bad_ind)) {
    stop(
      "independent_t_test rows must supply n1, n2 with t_df == n1 + n2 - 2 for: ",
      paste(wave1$claim_id[bad_ind], collapse = ", "),
      call. = FALSE
    )
  }

  os_rows <- wave1$frequentist_test %in%
    c("one_sample_t_test", "paired_t_test", "repeated_measures_anova")
  bad_os <- os_rows & !near_equal(wave1$n_total, wave1$df_for_bf + 1)
  if (any(bad_os)) {
    stop(
      "one-sample / paired / repeated-measures rows must have n_total == df + 1 ",
      "(use complete pairs, not consented N) for: ",
      paste(wave1$claim_id[bad_os], collapse = ", "),
      call. = FALSE
    )
  }

  families_needed <- unique(wave1$bf_family)
  if (any(priors$rscale <= 0, na.rm = TRUE) || any(is.na(priors$rscale))) {
    stop("config/priors_wave1.csv contains a missing or non-positive rscale.", call. = FALSE)
  }
  prior_dupes <- priors |>
    count(family, prior_label) |>
    filter(n > 1)
  if (nrow(prior_dupes) > 0) {
    stop("Duplicate (family, prior_label) rows in config/priors_wave1.csv: ",
         paste(prior_dupes$family, prior_dupes$prior_label, sep = "/", collapse = ", "),
         call. = FALSE)
  }
  for (fam in families_needed) {
    if (!any(priors$family == fam)) {
      stop("config/priors_wave1.csv has no priors for required family: ", fam, call. = FALSE)
    }
    n_primary <- sum(priors$family == fam & priors$prior_label == "primary")
    if (n_primary != 1) {
      stop("Family ", fam, " must have exactly one primary prior in config/priors_wave1.csv; found ",
           n_primary, ".", call. = FALSE)
    }
  }

  bf_for_t_row <- function(t_val, df_val, n1, n2, rscale, design, direction = NA_character_) {
    null_interval <- if (is.na(direction)) NULL else
      if (identical(direction, "positive")) c(0, Inf) else c(-Inf, 0)
    if (identical(design, "independent")) {
      bf <- meta.ttestBF(t = t_val, n1 = n1, n2 = n2, rscale = rscale,
                         nullInterval = null_interval)
    } else {
      bf <- meta.ttestBF(t = t_val, n1 = df_val + 1, rscale = rscale,
                         nullInterval = null_interval)
    }
    tab <- extractBF(bf)
    c(bf = as.numeric(tab$bf)[1], error = as.numeric(tab$error)[1])
  }

  bf_for_r_row <- function(r_val, n_val, rscale, direction = NA_character_) {
    if (is.na(r_val) || abs(r_val) >= 1) {
      stop("Correlation r must satisfy -1 < r < 1; received ", r_val, call. = FALSE)
    }
    if (is.na(n_val) || n_val < 4) {
      stop("Correlation requires n >= 4; received ", n_val, call. = FALSE)
    }
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(20260720L)
    x <- rnorm(n_val); x <- (x - mean(x)) / sd(x)
    z <- rnorm(n_val); z <- (z - mean(z)) / sd(z)
    z_orth <- residuals(lm(z ~ x)); z_orth <- (z_orth - mean(z_orth)) / sd(z_orth)
    y <- r_val * x + sqrt(1 - r_val^2) * z_orth
    achieved <- stats::cor(x, y)
    if (!isTRUE(all.equal(achieved, r_val, tolerance = 1e-8))) {
      stop("Synthetic correlation construction failed: requested r = ", r_val,
           " but achieved r = ", achieved, call. = FALSE)
    }
    null_interval <- if (is.na(direction)) NULL else
      if (identical(direction, "positive")) c(0, 1) else c(-1, 0)
    bf <- correlationBF(x, y, rscale = rscale, nullInterval = null_interval)
    tab <- extractBF(bf)
    c(bf = as.numeric(tab$bf)[1], error = as.numeric(tab$error)[1])
  }

  results <- pmap_dfr(wave1, function(...) {
    row <- list(...)
    map_dfr(seq_len(nrow(priors)), function(i) {
      p <- priors[i, ]
      if (p$family != row$bf_family) return(NULL)
      dir_for_bf <- if (identical(row$bf_sidedness, "one_sided")) {
        row$bf_direction
      } else {
        NA_character_
      }
      res <- if (row$bf_family == "correlation") {
        bf_for_r_row(row$r_value, row$n_total, p$rscale, direction = dir_for_bf)
      } else {
        bf_for_t_row(row$t_for_bf, row$df_for_bf, row$n1, row$n2, p$rscale,
                     row$design, direction = dir_for_bf)
      }
      tibble(
        claim_id = row$claim_id,
        study_id = row$study_id,
        stat_test = row$frequentist_test,
        bf_family = row$bf_family,
        design = row$design,
        bf_sidedness = row$bf_sidedness,
        bf_direction = row$bf_direction,
        observed_sign = row$observed_sign,
        direction_matches_observed = row$direction_matches_observed,
        prior_label = p$prior_label,
        rscale = p$rscale,
        t_for_bf = row$t_for_bf,
        df_for_bf = row$df_for_bf,
        r_value = row$r_value,
        n1 = row$n1,
        n2 = row$n2,
        n_total = row$n_total,
        p_value = row$p_value,
        bf10 = res[["bf"]],
        bf_error = res[["error"]]
      )
    })
  })

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  write_csv(results, output_path)
  invisible(results)
}
