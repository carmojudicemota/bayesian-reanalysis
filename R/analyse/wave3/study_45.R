study_45_claim_ids <- function() {"study_45_claim_01"}

load_study_45_data <- function(path = "data/raw/study_45/Untitled3.sav") {
  raw <- haven::read_sav(path)
  vars <- c("AVG_Aware", "AVG_Dev", "AVG_Explain", "AVG_Import")
  m <- as.matrix(raw[, vars])
  m <- m[stats::complete.cases(m), , drop = FALSE]
  m
}

study_45_rtnorm <- function(mean, sd, lower, upper) {
  lo <- stats::pnorm(lower, mean, sd)
  hi <- stats::pnorm(upper, mean, sd)
  if (hi - lo < 1e-12) return(min(max(mean, lower), upper))
  stats::qnorm(stats::runif(1, lo, hi), mean, sd)
}

study_45_contrasts <- function() {
  C <- stats::contr.helmert(4)
  apply(C, 2, function(col) col / sqrt(sum(col^2)))
}

study_45_gibbs <- function(m, n_samples = 4000, n_burnin = 1000) {
  n <- nrow(m)
  k <- ncol(m)
  z <- t(apply(m, 1, function(row) stats::qnorm(rank(row, ties.method = "average") / (k + 1))))
  mu <- colMeans(z)
  mu <- mu - mean(mu)
  Cmat <- study_45_contrasts()
  draws <- matrix(0, n_samples, k - 1)
  for (t in seq_len(n_burnin + n_samples)) {
    for (i in seq_len(n)) {
      row <- m[i, ]
      for (j in seq_len(k)) {
        smaller <- z[i, row < row[j]]
        larger <- z[i, row > row[j]]
        lower <- if (length(smaller)) max(smaller) else -Inf
        upper <- if (length(larger)) min(larger) else Inf
        z[i, j] <- study_45_rtnorm(mu[j], 1, lower, upper)
      }
    }
    mu <- stats::rnorm(k, colMeans(z), sqrt(1 / n))
    mu <- mu - mean(mu)
    if (t > n_burnin) draws[t - n_burnin, ] <- as.numeric(crossprod(Cmat, mu))
  }
  draws
}

study_45_afbf <- function(draws, n) {
  J <- ncol(draws)
  m_c <- colMeans(draws)
  S_c <- stats::cov(draws)
  b <- J / n
  log_post0 <- mvtnorm::dmvnorm(rep(0, J), mean = m_c, sigma = S_c, log = TRUE)
  log_prior0 <- mvtnorm::dmvnorm(rep(0, J), mean = rep(0, J), sigma = S_c / b, log = TRUE)
  exp(log_prior0 - log_post0)
}

compute_study_45_bayes_factors <- function(claim, priors = NULL, n_samples = 4000, n_burnin = 1000) {
  m <- load_study_45_data()
  set.seed(123)
  draws <- study_45_gibbs(m, n_samples, n_burnin)
  bf10 <- study_45_afbf(draws, nrow(m))
  wave3_row(
    claim = claim,
    bf10 = bf10,
    model_null = "equal dimension means (Friedman exchangeability)",
    model_alt = "at least one dimension mean differs",
    bf_family = "rank_latent_normal_rm",
    prior_family = "adjusted_fractional",
    method = "latent_normal_friedman_gibbs_afbf",
    prior_label = "primary"
  )
}