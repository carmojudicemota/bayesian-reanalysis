source("scripts/configure_runtime.R")
configure_runtime()

source("R/analyse/wave3/wave3_helpers.R")
source("R/analyse/wave3/study_40.R")
source("R/analyse/wave2/study_37_full_glmm.R")
source("R/analyse/wave2/run_study_37_full_glmm.R")
source("R/analyse/wave2/study_55.R")

message("Running robustness: Study 40")
run_study_40_prior_sweep()
study_40_pretrend_check()

message("Running robustness: Study 37")
run_study_37_gelman_cauchy_sensitivity()

message("Running robustness: Study 55")
study_55_welch_aafbf_sensitivity()

message("03_run_robustness complete.")
