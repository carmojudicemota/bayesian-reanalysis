source("R/analyse/wave3/wave3_helpers.R")
for (f in c("R/analyse/wave3/study_26.R", "R/analyse/wave3/study_39.R",
            "R/analyse/wave3/study_40.R", "R/analyse/wave3/study_45.R")) {
  if (file.exists(f)) source(f)
}

stability_wave_of <- function(study_id) {
  wave2 <- c("study_06", "study_10", "study_13", "study_29", "study_35",
             "study_37", "study_43", "study_47", "study_55", "study_60")
  wave3 <- c("study_03", "study_18", "study_20", "study_26", "study_39",
             "study_40", "study_44", "study_45")
  dplyr::case_when(study_id %in% wave2 ~ "Wave 2", study_id %in% wave3 ~ "Wave 3",
                   TRUE ~ NA_character_)
}


stability_source_of <- function(method) {
  dplyr::case_when(
    is.na(method) ~ NA_character_,
    grepl("BFpack|bain", method) ~ "analytic",
    grepl("bridge", method) ~ "bridge_sampling",
    grepl("RoBTT|blavaan", method) ~ "mcmc_convergence",
    grepl("brms_growth_curve_order_restricted_bf|brms_ordinal_order_restricted_bf", method) ~ "mcmc_convergence",
    grepl("gibbs", method) ~ "reseed_gibbs",
    grepl("bayesfactor|factorial|simple_effect", method) ~ "montecarlo_bayesfactor",
    TRUE ~ NA_character_
  )
}


stability_category <- function(x, threshold = log10(3)) {
  ifelse(x > threshold, "H1", ifelse(x < -threshold, "H0", "Inconclusive"))
}

stability_interval_distance <- function(lower, upper, threshold = log10(3)) {
  taus <- c(-threshold, threshold)
  mapply(function(lo, hi) {
    d <- vapply(taus, function(tau) {
      if (lo <= tau && tau <= hi) 0 else min(abs(lo - tau), abs(hi - tau))
    }, numeric(1))
    min(d)
  }, lower, upper)
}

stability_split_rhat <- function(chains) {
  n <- nrow(chains)
  half <- floor(n / 2)
  s <- cbind(chains[1:half, , drop = FALSE], chains[(n - half + 1):n, , drop = FALSE])
  m <- ncol(s)
  nn <- nrow(s)
  means <- colMeans(s)
  grand <- mean(means)
  b <- nn / (m - 1) * sum((means - grand)^2)
  w <- mean(apply(s, 2, stats::var))
  sqrt((((nn - 1) / nn) * w + b / nn) / w)
}

stability_reseed_samples <- function(label, bf_fun, seeds) {
  bf <- vapply(seeds, function(s) {
    set.seed(s)
    as.numeric(bf_fun())
  }, numeric(1))
  tibble::tibble(label = label, seed = seeds, bf10 = bf, log10_bf10 = log10(bf))
}

stability_reseed_mixing <- function(draw_fun, seeds) {
  chains <- vapply(seeds, function(s) {
    set.seed(s)
    draw_fun()
  }, numeric(length({
    set.seed(seeds[1])
    draw_fun()
  })))
  stability_split_rhat(chains)
}


reseed_study_26 <- function(seeds = 1:20, mixing_seeds = 101:104,
                            n_samples = 20000, n_burnin = 5000, r = 1 / sqrt(2)) {
  purrr::map_dfr(study_26_claim_ids(), function(cid) {
    differences <- load_study_26_data(cid)
    samples <- stability_reseed_samples(cid, function() {
      draws <- study_26_signed_rank_samples(differences, r, n_samples, n_burnin)
      study_26_savage_dickey(draws, r)
    }, seeds)
    rhat <- stability_reseed_mixing(function() {
      study_26_signed_rank_samples(differences, r, n_samples, n_burnin)
    }, mixing_seeds)
    dplyr::mutate(samples, study_id = "study_26", mcmc_rhat = rhat)
  })
}

reseed_study_39 <- function(seeds = 1:20, mixing_seeds = 101:104,
                            n_samples = 20000, n_burnin = 5000, kappa = 1) {
  prior0 <- study_39_prior_density_zero(kappa)
  purrr::map_dfr(study_39_claim_ids(), function(cid) {
    data <- load_study_39_data(cid)
    samples <- stability_reseed_samples(cid, function() {
      draws <- study_39_spearman_samples(data$x, data$y, n_samples, n_burnin, kappa)
      study_39_savage_dickey(draws, prior0)
    }, seeds)
    rhat <- stability_reseed_mixing(function() {
      study_39_spearman_samples(data$x, data$y, n_samples, n_burnin, kappa)
    }, mixing_seeds)
    dplyr::mutate(samples, study_id = "study_39", mcmc_rhat = rhat)
  })
}

summarise_reseed <- function(samples) {
  if (nrow(samples) == 0) return(tibble::tibble())
  samples |>
    dplyr::group_by(.data$study_id, .data$label) |>
    dplyr::summarise(
      source = "reseed_gibbs",
      n_seeds = dplyr::n(),
      central = mean(.data$log10_bf10),
      lower = min(.data$log10_bf10),
      upper = max(.data$log10_bf10),
      dispersion_log10 = stats::sd(.data$log10_bf10),
      mcmc_rhat = dplyr::first(.data$mcmc_rhat),
      ess_bulk_min = NA_real_,
      ess_tail_min = NA_real_,
      divergences = NA_real_,
      converged = TRUE,
      bridge_span_log10 = NA_real_,
      .groups = "drop"
    ) |>
    dplyr::rename(claim_id = "label")
}

harvest_wave3_brms_diagnostics <- function(study_id) {
  cache <- switch(
    study_id,
    study_40 = study_40_cache_path(),
    study_45 = study_45_cache_path(),
    stop("No direct brms diagnostic cache configured for ", study_id, ".", call. = FALSE)
  )
  
  if (!file.exists(cache)) {
    return(list(
      rhat_max = NA_real_, ess_bulk_min = NA_real_, ess_tail_min = NA_real_,
      divergences = NA_real_
    ))
  }
  
  tab <- utils::read.csv(cache, stringsAsFactors = FALSE)
  required <- c("rhat_max", "ess_bulk_min", "ess_tail_min", "divergences")
  if (!all(required %in% names(tab))) {
    return(list(
      rhat_max = NA_real_, ess_bulk_min = NA_real_, ess_tail_min = NA_real_,
      divergences = NA_real_
    ))
  }
  
  if ("prior_label" %in% names(tab) && any(tab$prior_label == "primary")) {
    tab <- tab[tab$prior_label == "primary", , drop = FALSE]
  }
  
  list(
    rhat_max = as.numeric(tab$rhat_max[[1]]),
    ess_bulk_min = as.numeric(tab$ess_bulk_min[[1]]),
    ess_tail_min = as.numeric(tab$ess_tail_min[[1]]),
    divergences = as.numeric(tab$divergences[[1]])
  )
}
harvest_bridge_stat <- function(study_id, claim_id, col) {
  candidates <- c(
    file.path("outputs/intermediate", paste0(study_id, "_bayes_factors.csv")),
    file.path("outputs/tables",paste0(study_id, "_bayes_factors.csv")),
    file.path("outputs/tables",paste0(study_id, "_full_glmm_bayes_factors.csv"))
  )
  cache <- candidates[file.exists(candidates)][1]
  if (is.na(cache)) return(NA_real_)
  tab <- utils::read.csv(cache, stringsAsFactors = FALSE)
  if (!col %in% names(tab)) return(NA_real_)
  row <- tab[tab$claim_id == claim_id, , drop = FALSE]
  if ("prior_label" %in% names(row) && any(row$prior_label == "primary", na.rm = TRUE)) {
    row <- row[row$prior_label == "primary", , drop = FALSE]
  }
  v <- suppressWarnings(as.numeric(row[[col]][1]))
  if (length(v) < 1) NA_real_ else v
}

harvest_pipeline_stability <- function(results_csv = "outputs/tables/bayes_factor_results.csv") {
  if (!file.exists(results_csv)) return(tibble::tibble())
  
  res <- readr::read_csv(results_csv, show_col_types = FALSE) |>
    dplyr::mutate(source = stability_source_of(.data$method)) |>
    dplyr::filter(!is.na(.data$source), .data$source != "reseed_gibbs")
  
  if ("prior_label" %in% names(res)) {
    res <- res |>
      dplyr::group_by(.data$claim_id) |>
      dplyr::mutate(
        keep = if (any(.data$prior_label == "primary", na.rm = TRUE))
          .data$prior_label == "primary" else dplyr::row_number() == 1L
      ) |>
      dplyr::ungroup() |>
      dplyr::filter(.data$keep)
  }
  
  diag40 <- harvest_wave3_brms_diagnostics("study_40")
  diag45 <- harvest_wave3_brms_diagnostics("study_45")
  
  res |>
    dplyr::rowwise() |>
    dplyr::mutate(
      central = .data$log10_bf10,
      bridge_sd_cache = if (.data$source == "bridge_sampling")
        harvest_bridge_stat(.data$study_id, .data$claim_id, "bridge_sd_log10") else NA_real_,
      bridge_span_log10 = if (.data$source == "bridge_sampling")
        harvest_bridge_stat(.data$study_id, .data$claim_id, "bridge_span_log10") else NA_real_,
      dispersion_log10 = dplyr::case_when(
        .data$source == "montecarlo_bayesfactor" ~ as.numeric(.data$bf_error) / log(10),
        .data$source == "bridge_sampling" ~ dplyr::coalesce(as.numeric(.data$bf_error), bridge_sd_cache),
        .data$source == "analytic" ~ 0,
        TRUE ~ NA_real_
      ),
      mcmc_rhat = dplyr::case_when(
        .data$study_id == "study_40" ~ diag40$rhat_max,
        .data$study_id == "study_45" ~ diag45$rhat_max,
        TRUE ~ NA_real_
      ),
      ess_bulk_min = dplyr::case_when(
        .data$study_id == "study_40" ~ diag40$ess_bulk_min,
        .data$study_id == "study_45" ~ diag45$ess_bulk_min,
        TRUE ~ NA_real_
      ),
      ess_tail_min = dplyr::case_when(
        .data$study_id == "study_40" ~ diag40$ess_tail_min,
        .data$study_id == "study_45" ~ diag45$ess_tail_min,
        TRUE ~ NA_real_
      ),
      divergences = dplyr::case_when(
        .data$study_id == "study_40" ~ diag40$divergences,
        .data$study_id == "study_45" ~ diag45$divergences,
        TRUE ~ NA_real_
      ),
      converged = dplyr::case_when(
        .data$study_id == "study_40" ~
          !is.na(diag40$rhat_max) && diag40$rhat_max < 1.01 &&
          !is.na(diag40$ess_bulk_min) && diag40$ess_bulk_min >= 400 &&
          !is.na(diag40$divergences) && diag40$divergences == 0,
        .data$study_id == "study_45" ~
          !is.na(diag45$rhat_max) && diag45$rhat_max < 1.01 &&
          !is.na(diag45$ess_bulk_min) && diag45$ess_bulk_min >= 400 &&
          !is.na(diag45$divergences) && diag45$divergences == 0,
        TRUE ~ TRUE
      ),
      lower = dplyr::if_else(
        is.na(.data$dispersion_log10), .data$central,
        .data$central - 2 * .data$dispersion_log10
      ),
      upper = dplyr::if_else(
        is.na(.data$dispersion_log10), .data$central,
        .data$central + 2 * .data$dispersion_log10
      ),
      n_seeds = NA_integer_
    ) |>
    dplyr::ungroup() |>
    dplyr::select(
      "claim_id", "study_id", "source", "n_seeds", "central", "lower", "upper",
      "dispersion_log10", "mcmc_rhat", "ess_bulk_min", "ess_tail_min",
      "divergences", "converged", "bridge_span_log10"
    )
}
build_stability_table <- function(reseed_summary, harvested,
                                  threshold = log10(3), magnitude_range_limit = 1,
                                  near_threshold_margin = 0.1) {
  combined <- dplyr::bind_rows(reseed_summary, harvested)
  if (nrow(combined) == 0) return(combined)
  combined |>
    dplyr::mutate(
      wave = stability_wave_of(.data$study_id),
      evidence_category = stability_category(.data$central, threshold),
      distance_to_threshold = pmin(abs(.data$central - threshold), abs(.data$central + threshold)),
      minimum_distance_to_threshold = stability_interval_distance(.data$lower, .data$upper, threshold),
      category_stable = stability_category(.data$lower, threshold) == stability_category(.data$upper, threshold),
      magnitude_reliable = dplyr::if_else(is.na(.data$dispersion_log10), NA, .data$dispersion_log10 < magnitude_range_limit),
      verdict = dplyr::case_when(
        .data$source == "analytic" ~ "Analytic (exact)",
        .data$source == "mcmc_convergence" & !.data$converged ~ "Unstable (no convergence)",
        !.data$category_stable ~ "Unstable",
        .data$minimum_distance_to_threshold < near_threshold_margin ~ "Borderline",
        !is.na(.data$magnitude_reliable) & !.data$magnitude_reliable ~ "Magnitude unreliable",
        TRUE ~ "Stable"
      ),
      stability_interpretation = dplyr::case_when(
        .data$source == "analytic" ~ "Closed-form Bayes factor; no Monte Carlo error",
        .data$verdict == "Unstable (no convergence)" ~ "Sampler did not meet the convergence threshold",
        .data$verdict == "Unstable" ~ "Evidence category changes across the numerical spread",
        .data$verdict == "Borderline" ~ "Category stable but the spread reaches near a threshold",
        .data$verdict == "Magnitude unreliable" ~ "Category stable; Bayes-factor magnitude not reliably estimated",
        TRUE ~ "Category and magnitude stable"
      )
    ) |>
    dplyr::arrange(.data$distance_to_threshold)
}


plot_stability_forest <- function(table, out_png = "outputs/figures/stability_forest.png", threshold = log10(3)) {
  tr <- function(x) asinh(x)
  jeff_bf10 <- c(1/100, 1/30, 1/10, 1/3, 1, 3, 10, 30, 100, 1e4, 1e10, 1e30, 1e90)
  jeff_log10 <- log10(jeff_bf10)
  jeff_labels <- c("1/100", "1/30", "1/10", "1/3", "1", "3", "10", "30", "100",
                   "10^4", "10^10", "10^30", "10^90")
  jeff_lines <- log10(c(1/100, 1/30, 1/10, 1/3, 3, 10, 30, 100))
  plot_df <- table |>
    dplyr::mutate(claim_id = factor(.data$claim_id, levels = rev(unique(.data$claim_id))))
  xr <- range(c(plot_df$lower, plot_df$upper, plot_df$central), na.rm = TRUE) + c(-0.5, 0.5)
  keep <- jeff_log10 >= xr[1] & jeff_log10 <= xr[2] & !(jeff_labels %in% c("1/30", "30"))
  vl <- jeff_lines[jeff_lines >= xr[1] & jeff_lines <= xr[2]]
  verdict_colours <- c(Stable = "#2166AC", `Analytic (exact)` = "#4393C3",
                       Borderline = "#E08214", `Magnitude unreliable` = "#B2182B",
                       Unstable = "#67001F", `Unstable (no convergence)` = "#67001F")
  plot <- ggplot2::ggplot(plot_df, ggplot2::aes(tr(.data$central), .data$claim_id,
                                                colour = .data$verdict)) +
    ggplot2::annotate("rect", xmin = tr(-threshold), xmax = tr(threshold), ymin = -Inf, ymax = Inf, fill = "#9E9E9E", alpha = 0.15) +
    ggplot2::geom_vline(xintercept = tr(vl), colour = "grey80", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = c(tr(-threshold), tr(threshold)), linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = tr(.data$lower), xmax = tr(.data$upper)), orientation = "y", width = 0.25, linewidth = 0.5) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = tr(jeff_log10[keep]), labels = jeff_labels[keep],
                                sec.axis = ggplot2::sec_axis(~ ., breaks = tr(jeff_log10[keep]),
                                                             labels = sprintf("%.2f", jeff_log10[keep]),
                                                             name = expression(log[10](BF[10])))) +
    ggplot2::scale_colour_manual(values = verdict_colours, name = "Stability verdict") +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$wave), scales = "free_y", space = "free_y") +
    ggplot2::labs(x = expression(BF[10] ~ "(Jeffreys scale; bars = numerical spread)"),
                  y = NULL, title = "Project-Wide Monte Carlo Stability of Bayes Factors",
                  subtitle = "Grey Band = inconclusive region; Bars = ±2 SD, Seed Range, or Bridge Repetitions") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
                   axis.text.x.bottom = ggplot2::element_text(angle = 45, hjust = 1, size = 8))
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_png, plot, width = 10, height = 2.5 + 0.4 * nrow(table), dpi = 150)
  invisible(plot)
}

plot_reseed_detail <- function(samples, table,
                               out_png = "outputs/figures/wave3_mc_stability.png",
                               threshold = log10(3)) {
  if (nrow(samples) == 0) return(invisible(NULL))
  verdicts <- dplyr::select(table, claim_id = "claim_id", "verdict")
  annotated <- dplyr::left_join(samples, verdicts, by = c("label" = "claim_id"))
  means <- annotated |>
    dplyr::group_by(.data$label) |>
    dplyr::summarise(mean_log10 = mean(.data$log10_bf10), .groups = "drop")
  verdict_colours <- c(Stable = "#2166AC", Borderline = "#E08214",
                       `Magnitude unreliable` = "#B2182B", Unstable = "#67001F")
  plot <- ggplot2::ggplot(annotated, ggplot2::aes(.data$log10_bf10, .data$label)) +
    ggplot2::annotate("rect", xmin = -threshold, xmax = threshold, ymin = -Inf, ymax = Inf,
                      fill = "#9E9E9E", alpha = 0.15) +
    ggplot2::geom_vline(xintercept = c(-threshold, threshold), linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_jitter(ggplot2::aes(colour = .data$verdict), height = 0.12, width = 0,
                         alpha = 0.6, size = 1.7) +
    ggplot2::geom_point(data = means, ggplot2::aes(.data$mean_log10, .data$label),
                        inherit.aes = FALSE, shape = 18, size = 3, colour = "black") +
    ggplot2::scale_colour_manual(values = verdict_colours, name = "Stability verdict") +
    ggplot2::labs(x = expression(log[10](BF[10]) ~ "across seeds"), y = NULL,
                  title = "Wave 3 Rank Based Bayes Factors: Per-Seed Reseeding",
                  subtitle = "Diamonds = Across-Seed Mean") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_png, plot, width = 8, height = 2.6 + 0.5 * dplyr::n_distinct(annotated$label),
                  dpi = 150)
  invisible(plot)
}

refresh_stability_forest <- function(
    table_csv = "outputs/diagnostics/stability_table.csv",
    results_csv = "outputs/tables/bayes_factor_results.csv",
    out_png = "outputs/figures/stability_forest.png") {
  
  old_reseed <- tibble::tibble()
  if (file.exists(table_csv)) {
    old <- readr::read_csv(table_csv, show_col_types = FALSE)
    old_reseed <- old |>
      dplyr::filter(.data$source == "reseed_gibbs", .data$study_id %in% c("study_26", "study_39")) |>
      dplyr::transmute(
        claim_id = .data$claim_id,
        study_id = .data$study_id,
        source = .data$source,
        n_seeds = .data$n_seeds,
        central = .data$central,
        lower = .data$lower,
        upper = .data$upper,
        dispersion_log10 = .data$dispersion_log10,
        mcmc_rhat = .data$mcmc_rhat,
        ess_bulk_min = NA_real_,
        ess_tail_min = NA_real_,
        divergences = NA_real_,
        converged = .data$converged,
        bridge_span_log10 = .data$bridge_span_log10
      )
  }
  
  harvested <- harvest_pipeline_stability(results_csv)
  table <- build_stability_table(old_reseed, harvested)
  
  dir.create(dirname(table_csv), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(table, table_csv)
  plot_stability_forest(table, out_png)
  
  message(
    "Refreshed stability table/forest from current canonical results and cached diagnostics; ",
    "retained stored Study 26/39 reseed summaries without rerunning them."
  )
  invisible(table)
}


run_project_stability <- function(results_csv = "outputs/tables/bayes_factor_results.csv",
                                  reseed_seeds = 1:20, mixing_seeds = 101:104,
                                  out_csv = "outputs/diagnostics/stability_table.csv",
                                  out_forest = "outputs/figures/stability_forest.png",
                                  out_reseed = "outputs/figures/wave3_mc_stability.png") {
  samples <- dplyr::bind_rows(
    reseed_study_26(reseed_seeds, mixing_seeds),
    reseed_study_39(reseed_seeds, mixing_seeds)
  )
  reseed_summary <- summarise_reseed(samples)
  harvested <- harvest_pipeline_stability(results_csv)
  table <- build_stability_table(reseed_summary, harvested)
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(table, out_csv)
  plot_stability_forest(table, out_forest)
  plot_reseed_detail(samples, table, out_reseed)
  message("Project stability: ", sum(table$verdict %in% c("Stable", "Analytic (exact)")),
          "/", nrow(table), " claims stable or exact; ",
          sum(table$verdict == "Borderline"), " borderline, ",
          sum(grepl("Unstable|unreliable", table$verdict)), " flagged.")
  list(samples = samples, table = table)
}


