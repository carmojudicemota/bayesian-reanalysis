study_39_claim_ids <- function() {
  c("study_39_claim_01","study_39_claim_02")
}


load_study_39_raw <- function(path = "data/raw/study_39/Open_Pedagogy_Student_Perceptions.sav") {
  if (!file.exists(path)) {path <- "data/raw/study_39/Untitled3.sav"}
  haven::read_sav(path)
}


study_39_find_columns <- function(raw) {
  labels <- vapply(raw,
    function(x) {
      label <- attr(x, "label")
      if (is.null(label)) {""} else {as.character(label)}
    },
    character(1)
  )
  
  find_column <- function(name_pattern = NULL,label_pattern = NULL) {
    matches <- character()
    
    if (!is.null(name_pattern)) {
      matches <- c(matches,
                   grep(name_pattern,
                        names(raw),
                        value = TRUE,
                        ignore.case = TRUE
                        )
                   )
    }
    
    if (!is.null(label_pattern)) {
      matches <- c(matches,
                   names(raw)[grepl(label_pattern,labels,ignore.case = TRUE)]
                   )
    }
    unique(matches)[[1]]
  }
  
  list(
    motivation = find_column(
      name_pattern = "^HOWMOTIV",
      label_pattern = "motivating.*final product.*openly available"
    ),
    diversity = find_column(
      name_pattern = "^IBELIEV",
      label_pattern = "photographs add diversity"
    ),
    engagement = find_column(
      name_pattern = "^DIDTHISC",
      label_pattern = "course seem more engaging.*project"
    ),
    flickr = find_column(
      name_pattern = "^(V41|V41_A)$",
      label_pattern = "sharing your photos.*world.*Flickr"
    )
  )
}


study_39_numeric <- function(x) {
  x <- haven::zap_labels(x)
  if (is.factor(x)) {x <- as.character(x)}
  suppressWarnings(as.numeric(x))
}


load_study_39_data <- function(claim_id,path = "data/raw/study_39/Open_Pedagogy_Student_Perceptions.sav") {
  raw <- load_study_39_raw(path)
  columns <- study_39_find_columns(raw)
  variables <- switch(claim_id,
                      study_39_claim_01 = c(columns$motivation,
                                            columns$diversity),
                      study_39_claim_02 = c(columns$engagement,
                                            columns$flickr),
    
    stop("Unknown Study 39 claim: ", claim_id, call. = FALSE)
  )
  
  x <- study_39_numeric(raw[[variables[[1]]]])
  y <- study_39_numeric(raw[[variables[[2]]]])
  keep <- stats::complete.cases(x, y)
  tibble::tibble(x = x[keep],y = y[keep])
}


prepare_study_39_rank_data <- function(data) {
  x_levels <- sort(unique(data$x))
  y_levels <- sort(unique(data$y))
  
  list(
    N = nrow(data),
    K_x = length(x_levels),
    K_y = length(y_levels),
    x_category = match(data$x, x_levels),
    y_category = match(data$y, y_levels),
    x_levels = x_levels,
    y_levels = y_levels
  )
}

study_39_rtnorm <- function(mean, sd, lower, upper) {
  lo <- stats::pnorm(lower, mean, sd)
  hi <- stats::pnorm(upper, mean, sd)
  if (hi - lo < 1e-12) return(min(max(mean, lower), upper))
  stats::qnorm(stats::runif(1, lo, hi), mean, sd)
}

study_39_update_latent <- function(values, z, z_other, rho) {
  sdev <- sqrt(1 - rho^2)
  for (i in seq_along(values)) {
    smaller <- z[values < values[i]]
    larger <- z[values > values[i]]
    lower <- if (length(smaller)) max(smaller) else -Inf
    upper <- if (length(larger)) min(larger) else Inf
    z[i] <- study_39_rtnorm(rho * z_other[i], sdev, lower, upper)
  }
  z
}

study_39_update_rho <- function(zx, zy, rho, kappa, step = 0.3) {
  n <- length(zx)
  sxy <- sum(zx * zy)
  sxx <- sum(zx^2)
  syy <- sum(zy^2)
  target <- function(r) {
    -n / 2 * log(1 - r^2) - (sxx - 2 * r * sxy + syy) / (2 * (1 - r^2)) +
      (1 / kappa - 1) * log(1 - r^2) + log(1 - r^2)
  }
  proposal <- tanh(atanh(rho) + stats::rnorm(1, 0, step))
  if (log(stats::runif(1)) < target(proposal) - target(rho)) proposal else rho
}

study_39_spearman_samples <- function(x, y, n_samples = 5000, n_burnin = 1000, kappa = 1) {
  n <- length(x)
  zx <- stats::qnorm(rank(x, ties.method = "average") / (n + 1))
  zy <- stats::qnorm(rank(y, ties.method = "average") / (n + 1))
  rho <- stats::cor(zx, zy)
  draws <- numeric(n_samples)
  for (t in seq_len(n_burnin + n_samples)) {
    zx <- study_39_update_latent(x, zx, zy, rho)
    zy <- study_39_update_latent(y, zy, zx, rho)
    rho <- study_39_update_rho(zx, zy, rho, kappa)
    if (t > n_burnin) draws[t - n_burnin] <- rho
  }
  draws
}

study_39_prior_density_zero <- function(kappa) {
  a <- 1 / kappa
  0.5 * (0.25)^(a - 1) / beta(a, a)
}

study_39_savage_dickey <- function(draws, prior_density_zero) {
  fit <- logspline::logspline(draws, lbound = -1, ubound = 1)
  prior_density_zero / logspline::dlogspline(0, fit)
}

compute_study_39_bayes_factors <- function(claim, priors = NULL, n_samples = 20000, n_burnin = 5000) {
  data <- load_study_39_data(claim$claim_id)
  kappas <- c(narrow = 0.5, primary = 1, wide = 2)
  purrr::map_dfr(names(kappas), function(lbl) {
    kappa <- kappas[[lbl]]
    set.seed(123)
    draws <- study_39_spearman_samples(data$x, data$y, n_samples, n_burnin, kappa)
    wave3_row(
      claim = claim,
      bf10 = study_39_savage_dickey(draws, study_39_prior_density_zero(kappa)),
      model_null = "latent rho = 0",
      model_alt = "latent rho != 0",
      bf_family = "rank_latent_normal",
      prior_family = "stretched_beta",
      method = "latent_normal_spearman_gibbs",
      prior_label = lbl
    )
  })
}

