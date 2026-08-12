source("scripts/configure_runtime.R")
configure_runtime()

message("Running reproductions...")
source("scripts/run_all_reproductions.R")

message("Preparing data tables...")
source("R/prepare/build_verified_results.R")
build_verified_results_draft()
source("R/prepare/build_claims.R")
build_claims_draft()
file.copy("data/derived/claims_draft.csv", "data/derived/claims.csv", overwrite = TRUE)
message("Data prepared: claims_draft.csv -> claims.csv")

message("Starting wave 2 run")


source("R/analyse/wave2/study_13.R")
study13 <- run_study_13_bayes_factors()
dir.create("outputs/intermediate", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(study13, "outputs/intermediate/study_13_bayes_factors.csv")

message("Running study 37: full glmm")
source("R/analyse/wave2/study_37_full_glmm.R")
source("R/analyse/wave2/run_study_37_full_glmm.R")
run_study_37_full_glmm()

message("Running study 43")
source("R/analyse/wave2/study_43.R")
run_study_43_bayes_factors(
  output_path = "outputs/intermediate/study_43_bayes_factors.csv"
)

message("Running study 55")
source("R/analyse/wave2/study_55.R")
unlink("outputs/intermediate/study_55_models", recursive = TRUE, force = TRUE)
study55 <- run_study_55_bayes_factors()
readr::write_csv(study55, "outputs/intermediate/study_55_bayes_factors.csv")


message("Starting wave 3 run")


source("scripts/04_fit_wave3_direct.R")

message("Running robustness... ")
source("scripts/03_run_robustness.R")

message("Compiling analysis...")
source("scripts/01_run_analysis.R")

message("Compiling stability analysis...")
source("R/analyse/stability.R")
run_project_stability()

message("Compiling figures...")
source("scripts/02_run_figures.R")

message("Full pipeline ran sucessfully (yay)")


source("scripts/run_all_reproductions.R")
source("R/prepare/build_verified_results.R")
build_verified_results_draft()
source("R/prepare/build_claims.R")
build_claims_draft()
file.copy("data/derived/claims_draft.csv", "data/derived/claims.csv", overwrite = TRUE)
message("Promoted claims_draft.csv -> claims.csv")

message("Starting wave 2 run")


source("R/analyse/wave2/study_13.R")
study13 <- run_study_13_bayes_factors()
dir.create("outputs/intermediate", recursive = TRUE, showWarnings = FALSE)
readr::write_csv(study13, "outputs/intermediate/study_13_bayes_factors.csv")

message("Running study 37: full glmm")
source("R/analyse/wave2/study_37_full_glmm.R")
source("R/analyse/wave2/run_study_37_full_glmm.R")
run_study_37_full_glmm()

message("Running study 43")
source("R/analyse/wave2/study_43.R")
run_study_43_bayes_factors(
  output_path = "outputs/intermediate/study_43_bayes_factors.csv"
)

message("Running study 55")
source("R/analyse/wave2/study_55.R")
unlink("outputs/intermediate/study_55_models", recursive = TRUE, force = TRUE)
study55 <- run_study_55_bayes_factors()
readr::write_csv(study55, "outputs/intermediate/study_55_bayes_factors.csv")


message("Starting wave 3 run")


source("scripts/04_fit_wave3_direct.R")

message("Running robustness... ")
source("scripts/03_run_robustness.R")

message("Compiling analysis...")
source("scripts/01_run_analysis.R")

message("Compiling stability analysis...")
source("R/analyse/stability.R")
run_project_stability()

message("Compiling figures...")
source("scripts/02_run_figures.R")

message("Full pipeline ran sucessfully (yay)")
