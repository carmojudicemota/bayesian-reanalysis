study_26_claim_ids <- function() {
  c("study_26_claim_01", "study_26_claim_02")
}

load_study_26_data <- function(claim_id, path = "data/raw/study_26/Final_Data_Set_Spring_2022.sav") {
  raw <- haven::read_sav(path, user_na = FALSE)
  as_num <- function(v) as.numeric(haven::zap_labels(haven::zap_missing(v)))
  pair <- switch(
    claim_id,
    study_26_claim_01 = list(first = as_num(raw$womenrec2), second = as_num(raw$womenrec1)),
    study_26_claim_02 = list(first = as_num(raw$womenrec1), second = as_num(raw$menrec1)),
    stop("Unknown Study 26 claim: ", claim_id, call. = FALSE)
  )
  keep <- stats::complete.cases(pair$first, pair$second)
  differences <- pair$first[keep] - pair$second[keep]
  differences[differences != 0]
}

study_26_rtnorm <- function(mean, sd, lower, upper) {
  lo <- stats::pnorm(lower, mean, sd)
  hi <- stats::pnorm(upper, mean, sd)
  if (hi - lo < 1e-12) return(min(max(mean, lower), upper))
  stats::qnorm(stats::runif(1, lo, hi), mean, sd)
}

study_26_update_latent <- function(abs_diff, sign_diff, z, delta) {
  absz <- abs(z)
  for (i in seq_along(abs_diff)) {
    smaller <- absz[abs_diff < abs_diff[i]]
    larger <- absz[abs_diff > abs_diff[i]]
    lo_mag <- if (length(smaller)) max(smaller) else 0
    hi_mag <- if (length(larger)) min(larger) else Inf
    if (sign_diff[i] > 0) {
      z[i] <- study_26_rtnorm(delta, 1, lo_mag, hi_mag)
    } else {
      z[i] <- study_26_rtnorm(delta, 1, -hi_mag, -lo_mag)
    }
    absz[i] <- abs(z[i])
  }
  z
}

study_26_signed_rank_samples <- function(differences, r, n_samples = 5000, n_burnin = 1000) {
  n <- length(differences)
  abs_diff <- abs(differences)
  sign_diff <- sign(differences)
  mag <- stats::qnorm((rank(abs_diff, ties.method = "average") / (n + 1) + 1) / 2)
  z <- sign_diff * mag
  delta <- 0
  g <- r^2
  draws <- numeric(n_samples)
  for (t in seq_len(n_burnin + n_samples)) {
    z <- study_26_update_latent(abs_diff, sign_diff, z, delta)
    prec <- n + 1 / g
    delta <- stats::rnorm(1, sum(z) / prec, sqrt(1 / prec))
    g <- 1 / stats::rgamma(1, shape = 1, rate = (r^2 + delta^2) / 2)
    if (t > n_burnin) draws[t - n_burnin] <- delta
  }
  draws
}

study_26_savage_dickey <- function(draws, r) {
  fit <- logspline::logspline(draws)
  stats::dcauchy(0, 0, r) / logspline::dlogspline(0, fit)
}

compute_study_26_bayes_factors <- function(claim, priors = NULL, n_samples = 5000, n_burnin = 1000) {
  differences <- load_study_26_data(claim$claim_id)
  scales <- c(narrow = 0.5, primary = 1 / sqrt(2), wide = 1)
  purrr::map_dfr(names(scales), function(label) {
    set.seed(123)
    draws <- study_26_signed_rank_samples(differences, scales[[label]], n_samples, n_burnin)
    wave3_row(
      claim = claim,
      bf10 = study_26_savage_dickey(draws, scales[[label]]),
      model_null = "delta = 0",
      model_alt = "delta != 0",
      bf_family = "rank_latent_normal",
      prior_family = "cauchy",
      method = "latent_normal_signed_rank_gibbs",
      prior_label = label
    )
  })
}


