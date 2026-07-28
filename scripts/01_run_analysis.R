source("R/analyse/compute_wave1_bayes_factors.R")
source("R/analyse/compute_wave2_bayes_factors.R")
source("R/analyse/compute_wave3_bayes_factors.R")
source("R/analyse/combine_bayes_factor_results.R")
source("R/analyse/dataset.R")
source("R/analyse/classify_concordance.R")

compute_wave1_bayes_factors()     
compute_wave2_bayes_factors()  
compute_wave3_bayes_factors()
combine_bayes_factor_results()    
build_dataset_outputs()
build_concordance_outputs()
message("Analysis done")



