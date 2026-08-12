study_45_claim_ids <- function() {
  "study_45_claim_01"
}


study_45_ortho_helmert <- function() {
  apply(stats::contr.helmert(4), 2, function(x) x / sqrt(sum(x^2)))
}


study_45_items_long <- function(path = "data/raw/study_45/Psy_Lit_Data.csv") {
  raw <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  attrs <- c("SciUnder", "IndRes", "ResSk", "Evid", "MulPer", "SocCul", "AppTh", "RealWr?")
  raw_dims <- c(Aware = "Aware", Develop = "Dev", Explain = "Explain", Import = "Import")
  grid <- expand.grid(attribute = attrs, rawdim = names(raw_dims), stringsAsFactors = FALSE)
  grid$col <- paste0(grid$attribute, "_", grid$rawdim)
  grid$dimension <- raw_dims[grid$rawdim]
  grid <- grid[grid$col %in% names(raw), , drop = FALSE]
  rows <- lapply(seq_len(nrow(grid)), function(j) {
    data.frame(
      student = seq_len(nrow(raw)),
      rating = suppressWarnings(as.integer(raw[[grid$col[j]]])),
      dimension = grid$dimension[j],
      attribute = grid$attribute[j],
      stringsAsFactors = FALSE
    )
  })
  long <- do.call(rbind, rows)
  long <- long[stats::complete.cases(long) & long$rating %in% 1:7, ]
  long$dimension <- factor(long$dimension, levels = c("Aware", "Dev", "Explain", "Import"))
  stats::contrasts(long$dimension) <- study_45_ortho_helmert()
  long$rating <- factor(long$rating, ordered = TRUE, levels = 1:7)
  long
}


fit_study_45_ordinal <- function(long = study_45_items_long(), b_scale = 1, seed = 123,
                                 chains = 4,
                                 cores = min(chains, getOption("bayesian_reanalysis.cores", getOption("mc.cores", 1L))),
                                 adapt_delta = 0.95) {
  brms::brm(
    rating ~ dimension + (1 | student) + (1 | attribute),
    data = long, family = brms::cumulative("probit"),
    prior = brms::set_prior(paste0("normal(0, ", b_scale, ")"), class = "b"),
    chains = chains, cores = cores, iter = 4000, warmup = 1000, seed = seed, refresh = 0,
    control = list(adapt_delta = adapt_delta)
  )
}


study_45_order_bf <- function(fit, b_scale = 1, n_prior = 2e5, seed = 1) {
  C <- study_45_ortho_helmert()
  dr <- as.data.frame(brms::as_draws_df(fit))
  B <- cbind(dr[["b_dimension1"]], dr[["b_dimension2"]], dr[["b_dimension3"]])
  D <- B %*% t(C)
  post_ok <- D[, 4] > D[, 1] & D[, 4] > D[, 2] & D[, 4] > D[, 3] & D[, 2] < D[, 1] & D[, 2] < D[, 3]
  set.seed(seed)
  P <- matrix(stats::rnorm(3 * n_prior, 0, b_scale), ncol = 3)
  PD <- P %*% t(C)
  prior_ok <- PD[, 4] > PD[, 1] & PD[, 4] > PD[, 2] & PD[, 4] > PD[, 3] & PD[, 2] < PD[, 1] & PD[, 2] < PD[, 3]
  p <- mean(post_ok)
  q <- mean(prior_ok)
  fl <- 1 / (2 * length(post_ok))
  pc <- min(max(p, fl), 1 - fl)
  qc <- min(max(q, fl), 1 - fl)
  list(
    post_prop = p, prior_prop = q,
    bf_vs_unconstrained = p / q,
    bf_vs_complement = (pc / (1 - pc)) / (qc / (1 - qc)),
    rhat_max = max(posterior::summarise_draws(brms::as_draws(fit), "rhat")$rhat, na.rm = TRUE)
  )
}


study_45_cache_path <- function() {
  "outputs/intermediate/study_45_bayes_factors.csv"
}


run_study_45_ordinal_fit <- function(b_scale = 1, seed = 123, cache_path = study_45_cache_path()) {
  fit <- fit_study_45_ordinal(b_scale = b_scale, seed = seed)
  bf <- study_45_order_bf(fit, b_scale = b_scale)
  row <- data.frame(
    claim_id = "study_45_claim_01", prior_label = "primary",
    bf10 = bf$bf_vs_complement, post_prop = bf$post_prop,
    prior_prop = bf$prior_prop, rhat_max = bf$rhat_max
  )
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(row, cache_path, row.names = FALSE)
  message(sprintf("Wrote Study 45 ordinal BF to %s (bf10 = %.4g)", cache_path, bf$bf_vs_complement))
  invisible(row)
}


compute_study_45_bayes_factors <- function(claim, cache_path = study_45_cache_path()) {
  if (!file.exists(cache_path)) run_study_45_ordinal_fit(cache_path = cache_path)
  cached <- utils::read.csv(cache_path, stringsAsFactors = FALSE)
  rows <- cached[cached$claim_id == claim$claim_id, , drop = FALSE]
  if (nrow(rows) < 1) {
    stop("No cached Study 45 result for ", claim$claim_id, ". Run run_study_45_ordinal_fit() first.",
         call. = FALSE)
  }
  wave3_row(
    claim = claim,
    bf10 = rows$bf10[1],
    model_null = "no career-highest / development-lowest ordering across dimensions",
    model_alt = "Import > all others & Dev < all others (career importance highest, development lowest)",
    bf_family = "ordinal_order_restricted",
    prior_family = "brms_cumulative_probit_symmetric_encompassing",
    method = "brms_ordinal_order_restricted_bf",
    prior_label = "primary"
  )
}
