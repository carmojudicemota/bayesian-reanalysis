study_40_claim_ids <- function() {
  "study_40_claim_01"
}


load_study_40_data <- function(path = "data/raw/study_40/Wood_and_Cross_2024_Final_Data_deidentified.csv") {
  data <- utils::read.csv(path, stringsAsFactors = FALSE)
  data$opt_in <- factor(data$opt_in, levels = c(0, 1))
  outcomes <- paste0("p_tot", 1:5)
  data[outcomes] <- lapply(data[outcomes], as.numeric)
  data
}


study_40_long <- function(data = load_study_40_data()) {
  occ <- 1:5
  rows <- lapply(seq_len(nrow(data)), function(i) {
    data.frame(
      student = i,
      group = as.integer(as.character(data$opt_in[i])),
      occasion = occ,
      score = as.numeric(data[i, paste0("p_tot", occ)])
    )
  })
  long <- do.call(rbind, rows)
  long <- long[stats::complete.cases(long), ]
  long$score_z <- as.numeric(scale(long$score))
  long$post4 <- as.integer(long$occasion == 4)
  long$post5 <- as.integer(long$occasion == 5)
  long$occ_f <- factor(long$occasion)
  long
}


study_40_priors <- function(focal_scale = 0.5) {
  c(
    brms::set_prior("normal(0, 5)", class = "b"),
    brms::set_prior(paste0("normal(0, ", focal_scale, ")"), class = "b", coef = "group:post4"),
    brms::set_prior(paste0("normal(0, ", focal_scale, ")"), class = "b", coef = "group:post5")
  )
}


fit_study_40_growth <- function(long = study_40_long(), focal_scale = 0.5, seed = 123, cores = 4) {
  brms::brm(
    score_z ~ occ_f + group + group:post4 + group:post5 + (1 | student),
    data = long, family = gaussian(), prior = study_40_priors(focal_scale),
    chains = 4, cores = cores, iter = 4000, warmup = 1000, seed = seed, refresh = 0
  )
}


study_40_order_bf <- function(fit, focal_scale = 0.5, n_prior = 2e5, seed = 1) {
  dr <- as.data.frame(brms::as_draws_df(fit))
  b4 <- dr[["b_group:post4"]]
  b5 <- dr[["b_group:post5"]]
  post_ok <- b4 > 0 & b5 > 0
  set.seed(seed)
  prior_ok <- stats::rnorm(n_prior, 0, focal_scale) > 0 & stats::rnorm(n_prior, 0, focal_scale) > 0
  p <- mean(post_ok)
  q <- mean(prior_ok)
  fl <- 1 / (2 * length(post_ok))
  pc <- min(max(p, fl), 1 - fl)
  qc <- min(max(q, fl), 1 - fl)
  list(
    post_prop = p, prior_prop = q,
    bf_vs_unconstrained = p / q,
    bf_vs_complement = (pc / (1 - pc)) / (qc / (1 - qc)),
    delta4 = mean(b4), delta5 = mean(b5),
    rhat_max = max(posterior::summarise_draws(brms::as_draws(fit), "rhat")$rhat, na.rm = TRUE)
  )
}


study_40_cache_path <- function() {
  "outputs/intermediate/study_40_bayes_factors.csv"
}


run_study_40_direct <- function(focal_scale = 0.5, seed = 123, cache_path = study_40_cache_path()) {
  fit <- fit_study_40_growth(focal_scale = focal_scale, seed = seed)
  bf <- study_40_order_bf(fit, focal_scale = focal_scale)
  row <- data.frame(
    claim_id = "study_40_claim_01", prior_label = "primary",
    bf10 = bf$bf_vs_complement, delta4 = bf$delta4, delta5 = bf$delta5,
    post_prop = bf$post_prop, rhat_max_m7 = bf$rhat_max, rhat_max_m8 = bf$rhat_max
  )
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(row, cache_path, row.names = FALSE)
  message(sprintf("Wrote Study 40 direct BF to %s (bf10 = %.4g)", cache_path, bf$bf_vs_complement))
  invisible(row)
}


compute_study_40_bayes_factors <- function(claim, cache_path = study_40_cache_path()) {
  if (!file.exists(cache_path)) run_study_40_direct(cache_path = cache_path)
  cached <- utils::read.csv(cache_path, stringsAsFactors = FALSE)
  rows <- cached[cached$claim_id == claim$claim_id, , drop = FALSE]
  if (nrow(rows) < 1) {
    stop("No cached Study 40 result for ", claim$claim_id, ". Run run_study_40_direct() first.",
         call. = FALSE)
  }
  wave3_row(
    claim = claim,
    bf10 = rows$bf10[1],
    model_null = "no additional post-intervention gain for opt-in (delta4 <= 0 or delta5 <= 0)",
    model_alt = "opt-in gains post-intervention (delta4 > 0 & delta5 > 0)",
    bf_family = "growth_order_restricted",
    prior_family = "brms_growth_normal_focal_encompassing",
    method = "brms_growth_curve_order_restricted_bf",
    prior_label = "primary"
  )
}


run_study_40_prior_sweep <- function(scales = c(narrow = 0.25, primary = 0.5, wide = 1), seed = 123) {
  long <- study_40_long()
  rows <- purrr::map_dfr(seq_along(scales), function(i) {
    s <- as.numeric(scales[i])
    fit <- fit_study_40_growth(long, focal_scale = s, seed = seed)
    bf <- study_40_order_bf(fit, focal_scale = s)
    tibble::tibble(
      prior_label = names(scales)[i], focal_scale = s,
      delta4 = bf$delta4, delta5 = bf$delta5,
      post_prob = bf$post_prop, bf_vs_complement = bf$bf_vs_complement
    )
  })
  dir.create("outputs/intermediate", recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(rows, "outputs/intermediate/study_40_prior_sweep.csv")
  print(rows)
  invisible(rows)
}


study_40_long_pretrend <- function() {
  long <- study_40_long()
  long$pre_c <- pmin(long$occasion, 3)
  long$pre_c <- long$pre_c - mean(long$pre_c)
  long
}


study_40_pretrend_check <- function(focal_scale = 0.5, seed = 123, cores = 4) {
  long <- study_40_long_pretrend()
  pr <- c(
    brms::set_prior("normal(0, 5)", class = "b"),
    brms::set_prior(paste0("normal(0, ", focal_scale, ")"), class = "b", coef = "group:post4"),
    brms::set_prior(paste0("normal(0, ", focal_scale, ")"), class = "b", coef = "group:post5")
  )
  fit <- brms::brm(
    score_z ~ occ_f + group + group:pre_c + group:post4 + group:post5 + (1 | student),
    data = long, family = gaussian(), prior = pr,
    chains = 4, cores = cores, iter = 4000, warmup = 1000, seed = seed, refresh = 0
  )
  bf <- study_40_order_bf(fit, focal_scale = focal_scale)
  print(bf)
  invisible(bf)
}
