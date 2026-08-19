source("scripts/configure_runtime.R")
configure_runtime()

source("scripts/run_all_reproductions.R")
source("R/prepare/build_verified_results.R")
build_verified_results_draft()
source("R/prepare/build_claims.R")
build_claims_draft()
file.copy("data/derived/claims_draft.csv", "data/derived/claims.csv", overwrite = TRUE)
message("Promoted claims_draft.csv -> claims.csv")

source("R/analyse/wave2/study_37_full_glmm.R")
source("R/analyse/wave2/run_study_37_full_glmm.R")
run_study_37_full_glmm()

source("R/analyse/wave3/wave3_helpers.R")
source("R/analyse/wave3/study_40.R")
run_study_40_direct()
source("R/analyse/wave3/study_45.R")
run_study_45_ordinal_fit()

source("scripts/01_run_analysis.R")
source("scripts/02_run_figures.R")
message("00_run_all complete (all models refit from scratch).")
