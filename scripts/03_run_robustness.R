source("R/analyse/wave3/wave3_helpers.R")
source("R/analyse/wave3/study_40.R")
source("R/analyse/wave2/study_37_full_glmm.R")
source("R/analyse/wave2/run_study_37_full_glmm.R")
source("R/analyse/wave2/study_55.R")

message("Robustness suite: Study 40 loading sweep + WAIC/LOO ...")
run_study_40_bayes_factors()
study_40_predictive_comparison()
study_40_predictive_sweep()

message("Robustness suite: Study 37 Cauchy(0, 2.5) focal-prior check ...")
run_study_37_gelman_cauchy_sensitivity()

message("Robustness suite: Study 55 adjusted-fractional Bayesian Welch ANOVA ...")
study_55_welch_aafbf_sensitivity()

message("03_run_robustness complete. (MANOVA Student-t robustness runs within the main analysis.)")
