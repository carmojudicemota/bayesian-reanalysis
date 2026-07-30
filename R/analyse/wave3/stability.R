source("R/analyse/wave3/wave3_helpers.R")
for (f in c("R/analyse/wave3/study_26.R", "R/analyse/wave3/study_39.R", "R/analyse/wave3/study_45.R")) {
  if (file.exists(f)) source(f)
}


wave3_stability_samples <- function(label, bf_fun, seeds) {
  bf <- vapply(seeds, function(s) {
    set.seed(s)
    as.numeric(bf_fun())
  }, numeric(1))
  tibble::tibble(label = label, seed = seeds, bf10 = bf, log10_bf10 = log10(bf))
}


wave3_stability_category <- function(log10_bf10, threshold = log10(3)) {
  ifelse(log10_bf10 > threshold, "H1", ifelse(log10_bf10 < -threshold, "H0", "Inconclusive"))
}


wave3_stability_summary <- function(samples, threshold = log10(3)) {
  samples |>
    dplyr::group_by(.data$label) |>
    dplyr::summarise(
      n_seeds = dplyr::n(),
      mean_log10 = mean(.data$log10_bf10),
      sd_log10 = stats::sd(.data$log10_bf10),
      min_log10 = min(.data$log10_bf10),
      max_log10 = max(.data$log10_bf10),
      range_log10 = max(.data$log10_bf10) - min(.data$log10_bf10),
      distance_to_threshold = min(abs(mean(.data$log10_bf10) - c(-threshold, threshold))),
      category_stable = dplyr::n_distinct(wave3_stability_category(.data$log10_bf10, threshold)) == 1L,
      .groups = "drop"
    )
}


stability_study_26 <- function(seeds = 1:20, n_samples = 5000, n_burnin = 1000, r = 1 / sqrt(2)) {
  purrr::map_dfr(study_26_claim_ids(), function(cid) {
    differences <- load_study_26_data(cid)
    wave3_stability_samples(cid, function() {
      draws <- study_26_signed_rank_samples(differences, r, n_samples, n_burnin)
      study_26_savage_dickey(draws, r)
    }, seeds)
  })
}


stability_study_39 <- function(seeds = 1:20, n_samples = 5000, n_burnin = 1000, kappa = 1) {
  prior0 <- study_39_prior_density_zero(kappa)
  purrr::map_dfr(study_39_claim_ids(), function(cid) {
    data <- load_study_39_data(cid)
    wave3_stability_samples(cid, function() {
      draws <- study_39_spearman_samples(data$x, data$y, n_samples, n_burnin, kappa)
      study_39_savage_dickey(draws, prior0)
    }, seeds)
  })
}


stability_study_45 <- function(seeds = 1:8, n_samples = 4000, n_burnin = 1000) {
  if (!exists("study_45_gibbs")) {
    stop("Study 45 functions not found; source R/analyse/wave3/study_45.R first.", call. = FALSE)
  }
  m <- load_study_45_data()
  wave3_stability_samples("study_45_claim_01", function() {
    draws <- study_45_gibbs(m, n_samples, n_burnin)
    study_45_afbf(draws, nrow(m))
  }, seeds)
}


plot_wave3_stability <- function(samples, out_png = "outputs/figures/wave3_mc_stability.png",
                                 threshold = log10(3)) {
  summary <- wave3_stability_summary(samples, threshold)
  annotated <- dplyr::left_join(samples, dplyr::select(summary, "label", "category_stable"),
                                by = "label")
  plot <- ggplot2::ggplot(annotated, ggplot2::aes(.data$log10_bf10, .data$label)) +
    ggplot2::annotate("rect", xmin = -threshold, xmax = threshold, ymin = -Inf, ymax = Inf,
                      fill = "#9E9E9E", alpha = 0.15) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = c(-threshold, threshold), linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_jitter(ggplot2::aes(colour = .data$category_stable), height = 0.12, width = 0,
                         alpha = 0.6, size = 1.7) +
    ggplot2::geom_point(data = summary, ggplot2::aes(.data$mean_log10, .data$label),
                        inherit.aes = FALSE, shape = 18, size = 3, colour = "black") +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#B2182B"),
                                 name = "Evidence category stable across seeds") +
    ggplot2::labs(x = expression(log[10](BF[10]) ~ "across seeds"), y = NULL,
                  title = "Wave 3 Monte Carlo stability of rank-based Bayes factors",
                  subtitle = "Grey band = inconclusive region; diamonds = across-seed mean") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_png, plot, width = 8,
                  height = 2.6 + 0.5 * dplyr::n_distinct(annotated$label), dpi = 150)
  invisible(plot)
}


run_wave3_stability <- function(seeds = 1:20, seeds_45 = 1:8,
                                out_csv = "outputs/diagnostics/wave3_mc_stability.csv",
                                out_png = "outputs/figures/wave3_mc_stability.png") {
  parts <- list(stability_study_26(seeds), stability_study_39(seeds))
  if (exists("study_45_gibbs")) {
    parts <- c(parts, list(stability_study_45(seeds_45)))
  }
  samples <- dplyr::bind_rows(parts)
  summary <- wave3_stability_summary(samples)
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(summary, out_csv)
  plot_wave3_stability(samples, out_png)
  message("Wave 3 MC stability: ", sum(summary$category_stable), "/", nrow(summary),
          " claims category-stable across seeds.")
  list(samples = samples, summary = summary)
}
