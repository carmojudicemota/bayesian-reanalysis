configure_runtime <- function(
    max_cores = 4L,
    force = FALSE,
    verbose = TRUE) {
  
  if (!force) {
    current <- getOption("bayesian_reanalysis.cores", NULL)
    if (!is.null(current)) {
      return(invisible(as.integer(current)))
      }
  }
  
  detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (
    length(detected) != 1L ||
    is.na(detected) ||
    !is.finite(detected) ||
    detected < 1L
  ) {detected <- 1L}
  
  detected <- as.integer(detected)
  
  
  env_cores <- trimws(Sys.getenv("BAYES_REANALYSIS_CORES", unset = ""))
  if (nzchar(env_cores)) {
    
    requested <- suppressWarnings(as.integer(env_cores))
    
    if (is.na(requested) || requested < 1L) {
      stop("BAYES_REANALYSIS_CORES must be a positive integer; got: ",env_cores,call. = FALSE)
    }
    
  } else {requested <- as.integer(max_cores)}
  
  if (
    length(max_cores) != 1L ||
    is.na(max_cores) ||
    !is.finite(max_cores) ||
    max_cores < 1L
  ) {
    stop("max_cores must be a positive integer.",call. = FALSE)
  }
  
  max_cores <- as.integer(max_cores)
  
  n_cores <- max(
    1L,
    min(
      detected,
      requested,
      max_cores
    )
  )
  
  options(
    bayesian_reanalysis.cores = n_cores,
    mc.cores = n_cores
  )
  
  if (isTRUE(verbose)) {
    message(
      "Runtime configured: using up to ",
      n_cores,
      " parallel core",
      if (n_cores == 1L) "" else "s",
      " for MCMC/bridge computations (detected ",
      detected,
      ")."
    )
  }
  
  invisible(n_cores)
}