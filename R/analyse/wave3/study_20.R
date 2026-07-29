study_20_claim_ids <- function() {
  c("study_20_claim_01", "study_20_claim_02")
}


load_study_20_data <- function(path = "data/raw/study_20/PA_sj-csv-2-plj-10.1177_14757257221098860.csv") {
  dat <- utils::read.csv(path, sep = ";", header = TRUE, na.strings = c("-999", "", "NA"), stringsAsFactors = FALSE)
  mean_of <- function(cols) rowMeans(dat[, cols])
  dat$FRG_A_t1 <- mean_of(paste0("FRG_A_0", 1:5, "_t1"))
  dat$FRG_M_t1 <- mean_of(paste0("FRG_M_0", 1:5, "_t1"))
  dat$FRG_E_t1 <- mean_of(paste0("FRG_E_0", 1:5, "_t1"))
  dat$FRG_A_t2 <- mean_of(paste0("FRG_A_0", 1:5, "_t2"))
  dat$FRG_M_t2 <- mean_of(paste0("FRG_M_0", 1:5, "_t2"))
  dat$FRG_E_t2 <- mean_of(paste0("FRG_E_0", 1:5, "_t2"))
  dat$D_FRG_Index_t1 <- dat$FRG_E_t1 - 0.5 * (dat$FRG_A_t1 + dat$FRG_M_t1)
  dat$D_FRG_Index_t2 <- dat$FRG_E_t2 - 0.5 * (dat$FRG_A_t2 + dat$FRG_M_t2)
  pre_mean <- mean(dat$D_FRG_Index_t1, na.rm = TRUE)
  pre_sd <- stats::sd(dat$D_FRG_Index_t1, na.rm = TRUE)
  dat$D_FRG_Index_t2 <- (dat$D_FRG_Index_t2 - pre_mean) / pre_sd
  dat$D_FRG_Index_t1 <- (dat$D_FRG_Index_t1 - pre_mean) / pre_sd
  dat$DV <- as.numeric(dat$group_t2 != "Kontrolle")
  dat$R <- as.numeric(!(dat$group_t2 %in% c("Kontrolle", "Volition")))
  dat$Social_Interaction <- as.numeric(dat$group_t2 == "Soziale Interaktion")
  dat
}


study_20_model <- function() {
  '
    delta_FRG =~ 1*D_FRG_Index_t2
    D_FRG_Index_t2 ~ 1*D_FRG_Index_t1
    delta_FRG ~ b0*1
    D_FRG_Index_t1 ~ 1
    D_FRG_Index_t2 ~ 0
    delta_FRG ~ D_FRG_Index_t1
    delta_FRG ~~ delta_FRG
    D_FRG_Index_t1 ~~ D_FRG_Index_t1
    D_FRG_Index_t2 ~~ 0*D_FRG_Index_t2
    delta_FRG ~ b1*DV + b2*R + b3*Social_Interaction
    D_FRG_Index_t1 ~ DV + R + Social_Interaction
    m_C := b0
    m_DV := b0 + b1
    m_R := b0 + b1 + b2
    m_SI := b0 + b1 + b2 + b3
    Test_H1a := (m_DV + m_R + m_SI)/3 - b0
    Test_H1c := m_SI - m_R
  '
}


fit_study_20_model <- function(data) {
  lavaan::sem(model = study_20_model(), data = data, estimator = "ML", std.lv = FALSE, fixed.x = FALSE, missing = "FIML")
}


study_20_parameter <- function(claim_id) {
  switch(
    claim_id,
    study_20_claim_01 = "Test_H1a",
    study_20_claim_02 = "Test_H1c",
    stop("Unknown Study 20 claim: ", claim_id, call. = FALSE)
  )
}


study_20_standardized_contrast <- function(fit, par_name) {
  std <- lavaan::standardizedSolution(fit)
  row <- std[std$op == ":=" & std$lhs == par_name, ]
  if (nrow(row) != 1) {
    stop("Could not extract standardized contrast ", par_name, " from lavaan fit.", call. = FALSE)
  }
  list(estimate = as.numeric(row$est.std), se = as.numeric(row$se))
}


study_20_fraction_bf <- function(estimate, se, n, fraction) {
  theta <- c(psi = estimate)
  result <- bain::bain(theta, hypothesis = "psi = 0; psi > 0", n = n, Sigma = list(matrix(se^2, 1, 1)), group_parameters = 1, joint_parameters = 0, fraction = fraction)
  as.numeric(result$BFmatrix["H2", "H1"])
}


compute_study_20_bayes_factors <- function(claim, priors = NULL) {
  data <- load_study_20_data()
  fit <- fit_study_20_model(data)
  par_name <- study_20_parameter(claim$claim_id)
  contrast <- study_20_standardized_contrast(fit, par_name)
  n <- nrow(data)
  fractions <- c(primary = 1, sensitivity_b2 = 2, sensitivity_b3 = 3)
  set.seed(123)
  purrr::map_dfr(names(fractions), function(label) {
    bf10 <- study_20_fraction_bf(contrast$estimate, contrast$se, n, fractions[[label]])
    wave3_row(
      claim = claim,
      bf10 = bf10,
      model_null = paste0(par_name, " = 0"),
      model_alt = paste0(par_name, " > 0"),
      bf_family = "adjusted_fractional_order_restricted",
      prior_family = "fractional",
      method = "bain_aafbf_standardized_lavaan_contrast",
      prior_label = label
    )
  })
}
