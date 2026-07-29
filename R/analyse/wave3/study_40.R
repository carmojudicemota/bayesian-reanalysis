study_40_claim_ids <- function() {
  "study_40_claim_01"
}


load_study_40_data <- function(path = "data/raw/study_40/Wood_and_Cross_2024_Final_Data_deidentified.csv") {
  data <- utils::read.csv(path, stringsAsFactors = FALSE)
  data$opt_in <- factor(data$opt_in, levels = c(0, 1))
  outcomes <- paste0("p_tot", 1:5)
  data[outcomes] <- lapply(data[outcomes], as.numeric)
  scale_constant <- stats::sd(as.matrix(data[outcomes]), na.rm = TRUE)
  data[outcomes] <- data[outcomes] / scale_constant
  attr(data, "scale_constant") <- scale_constant
  data
}


study_40_model_7 <- function() {
  '
    i =~ 1*p_tot1 + 1*p_tot2 + 1*p_tot3 + 1*p_tot4 + 1*p_tot5
    s =~ 0*p_tot1 + m*p_tot2 + 1*p_tot3 + NA*p_tot4 + NA*p_tot5
  '
}


study_40_model_8 <- function() {
  '
    i =~ 1*p_tot1 + 1*p_tot2 + 1*p_tot3 + 1*p_tot4 + 1*p_tot5
    s =~ 0*p_tot1 + NA*p_tot2 + 1*p_tot3 + NA*p_tot4 + NA*p_tot5
  '
}


study_40_priors <- function(loading_scale = 0.5) {
  blavaan::dpriors(
    lambda = paste0("normal(1,", loading_scale, ")"),
    rho = "beta(2,2)"
  )
}


fit_study_40_model_7 <- function(data, burnin = 5000, sample = 10000, seed = 123,
                                 loading_scale = 0.5) {
  blavaan::bgrowth(
    model = study_40_model_7(),
    data = data,
    group = "opt_in",
    group.equal = c("residuals", "lv.variances", "lv.covariances"),
    cp = "srs",
    dp = study_40_priors(loading_scale),
    n.chains = 4,
    burnin = burnin,
    sample = sample,
    adapt = 2000,
    seed = seed,
    bcontrol = list(control = list(adapt_delta = 0.999, max_treedepth = 15))
  )
}


fit_study_40_model_8 <- function(data, burnin = 5000, sample = 10000, seed = 124,
                                 loading_scale = 0.5) {
  blavaan::bgrowth(
    model = study_40_model_8(),
    data = data,
    group = "opt_in",
    group.equal = c("loadings", "residuals", "lv.variances", "lv.covariances"),
    cp = "srs",
    dp = study_40_priors(loading_scale),
    n.chains = 4,
    burnin = burnin,
    sample = sample,
    adapt = 2000,
    seed = seed,
    bcontrol = list(control = list(adapt_delta = 0.999, max_treedepth = 15))
  )
}


study_40_diagnostics <- function(fit) {
  rhat <- blavaan::blavInspect(fit, "rhat")
  ess <- blavaan::blavInspect(fit, "neff")
  list(
    rhat_max = max(rhat, na.rm = TRUE),
    ess_bulk_min = min(ess, na.rm = TRUE),
    converged = max(rhat, na.rm = TRUE) < 1.01 && min(ess, na.rm = TRUE) > 400
  )
}
