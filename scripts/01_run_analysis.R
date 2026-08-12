source("scripts/configure_runtime.R")
configure_runtime()

source("R/analyse/compute_wave1_bayes_factors.R")
source("R/analyse/compute_wave2_bayes_factors.R")
source("R/analyse/compute_wave3_bayes_factors.R")
source("R/analyse/combine_bayes_factor_results.R")
source("R/analyse/classify_concordance.R")

message("Computing wave 1 Bayes factors")
compute_wave1_bayes_factors()
message("Computing wave 2 Bayes factors")
compute_wave2_bayes_factors()
message("Computing wave 3 Bayes factors")
compute_wave3_bayes_factors()
message("Combining Bayes factors results")
combine_bayes_factor_results()
message("Building concordance outputs")
build_concordance_outputs()

message("Analysing disagreement")
source("R/analyse/disagreement_analysis.R")

message("Analysis done")
