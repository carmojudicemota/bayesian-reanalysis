#runs everything from saved outputs, does not refit models
source("scripts/configure_runtime.R")
configure_runtime()

message("Compiling analysis...")
source("scripts/01_run_analysis.R")

message("Compiling figures...")
source("scripts/02_run_figures.R")


message("Pipeline ran sucessfully")
