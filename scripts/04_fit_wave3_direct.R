source("R/analyse/wave3/wave3_helpers.R")
source("R/analyse/wave3/study_40.R")
source("R/analyse/wave3/study_45.R")

message("Fitting Study 40 direct growth model (from scratch) ...")
run_study_40_direct()

message("Fitting Study 45 ordinal model (from scratch) ...")
run_study_45_ordinal_fit()

message("04_fit_wave3_direct complete.")
