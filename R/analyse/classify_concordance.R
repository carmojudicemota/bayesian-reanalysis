library(readr)
library(dplyr)
library(ggplot2)
library(forcats)

p_band_levels <- c("p >= .05", ".01 <= p < .05", ".001 <= p < .01", "p < .001")

bf_strength_levels <- c(
  "Decisive/extreme support for H0", "Very strong support for H0", "Strong support for H0",
  "Substantial/moderate support for H0", "Weak/anecdotal support for H0", "Equal support",
  "Weak/anecdotal support for H1", "Substantial/moderate support for H1", "Strong support for H1",
  "Very strong support for H1", "Decisive/extreme support for H1"
)

concordance_cell_levels <- c(
  "Significant + H1", "Significant + inconclusive", "Significant + H0",
  "Nonsignificant + H0", "Nonsignificant + inconclusive", "Nonsignificant + H1"
)

classify_p_band <- function(p) {
  result <- case_when(
    is.na(p) ~ NA_character_, p < .001 ~ "p < .001", p < .01 ~ ".001 <= p < .01",
    p < .05 ~ ".01 <= p < .05", TRUE ~ "p >= .05"
  )
  factor(result, levels = p_band_levels, ordered = TRUE)
}

classify_frequentist_result <- function(p, alpha = .05) {
  result <- case_when(is.na(p) ~ NA_character_, p < alpha ~ "Significant", TRUE ~ "Nonsignificant")
  factor(result, levels = c("Nonsignificant", "Significant"))
}

classify_bf_conclusion <- function(bf10, k = 3) {
  result <- case_when(is.na(bf10) ~ NA_character_, bf10 >= k ~ "H1", bf10 <= 1 / k ~ "H0", TRUE ~ "Inconclusive")
  factor(result, levels = c("H0", "Inconclusive", "H1"))
}

classify_bf_strength <- function(bf10) {
  result <- case_when(
    is.na(bf10) ~ NA_character_,
    # Boundaries are inclusive on the H0 side so the categories are exact
    # reciprocals of the H1 side: the H1 side reads [1,3) [3,10) [10,30) [30,100)
    # [100,Inf), so the H0 side must read (1/3,1) (1/10,1/3] (1/30,1/10]
    # (1/100,1/30] (0,1/100]. With "<" here, BF = 1/10 was labelled "Substantial"
    # while BF = 10 was labelled "Strong" -- the same evidence, different names.
    bf10 <= 1 / 100 ~ "Decisive/extreme support for H0",
    bf10 <= 1 / 30 ~ "Very strong support for H0",
    bf10 <= 1 / 10 ~ "Strong support for H0",
    bf10 <= 1 / 3 ~ "Substantial/moderate support for H0",
    bf10 < 1 ~ "Weak/anecdotal support for H0",
    bf10 == 1 ~ "Equal support",
    bf10 < 3 ~ "Weak/anecdotal support for H1",
    bf10 < 10 ~ "Substantial/moderate support for H1",
    bf10 < 30 ~ "Strong support for H1",
    bf10 < 100 ~ "Very strong support for H1",
    TRUE ~ "Decisive/extreme support for H1"
  )
  factor(result, levels = bf_strength_levels, ordered = TRUE)
}

classify_favoured_side <- function(bf10) {
  case_when(is.na(bf10) ~ NA_character_, bf10 > 1 ~ "H1", bf10 < 1 ~ "H0", TRUE ~ "Equal")
}

classify_concordance_cell <- function(frequentist_result, bf_conclusion) {
  result <- case_when(
    frequentist_result == "Significant" & bf_conclusion == "H1" ~ "Significant + H1",
    frequentist_result == "Significant" & bf_conclusion == "Inconclusive" ~ "Significant + inconclusive",
    frequentist_result == "Significant" & bf_conclusion == "H0" ~ "Significant + H0",
    frequentist_result == "Nonsignificant" & bf_conclusion == "H0" ~ "Nonsignificant + H0",
    frequentist_result == "Nonsignificant" & bf_conclusion == "Inconclusive" ~ "Nonsignificant + inconclusive",
    frequentist_result == "Nonsignificant" & bf_conclusion == "H1" ~ "Nonsignificant + H1",
    TRUE ~ NA_character_
  )
  factor(result, levels = concordance_cell_levels)
}

classify_concordance_status <- function(cell) {
  result <- case_when(
    cell %in% c("Significant + H1", "Nonsignificant + H0") ~ "Concordant",
    cell %in% c("Significant + inconclusive", "Nonsignificant + inconclusive") ~ "Inconclusive",
    cell %in% c("Significant + H0", "Nonsignificant + H1") ~ "Discordant",
    TRUE ~ NA_character_
  )
  factor(result, levels = c("Concordant", "Inconclusive", "Discordant"))
}

add_evidence_classifications <- function(results, alpha, k) {
  results |>
    mutate(
      p_band = classify_p_band(p_value),
      frequentist_result = classify_frequentist_result(p_value, alpha = alpha),
      bf_conclusion = classify_bf_conclusion(bf10, k = k),
      bf_strength = classify_bf_strength(bf10),
      favoured_side = classify_favoured_side(bf10),
      concordance_cell = classify_concordance_cell(frequentist_result, bf_conclusion),
      concordance_status = classify_concordance_status(concordance_cell),
      negative_log10_p = -log10(p_value),
      log10_bf10 = log10(bf10)
    )
}

summarise_prior_sensitivity <- function(results) {
  results |>
    group_by(claim_id, study_id) |>
    summarise(
      n_prior_specifications = n(),
      prior_labels = paste(sort(unique(prior_label)), collapse = "; "),
      minimum_bf10 = min(bf10),
      maximum_bf10 = max(bf10),
      minimum_log10_bf10 = min(log10_bf10),
      maximum_log10_bf10 = max(log10_bf10),
      prior_sensitivity_span = maximum_log10_bf10 - minimum_log10_bf10,
      bf_conclusion_changed = n_distinct(as.character(bf_conclusion)) > 1,
      favoured_side_changed = n_distinct(favoured_side) > 1,
      .groups = "drop"
    )
}

build_concordance_summary <- function(claim_level) {
  total_claims <- nrow(claim_level)
  
  claim_counts <- claim_level |>
    mutate(concordance_cell = as.character(concordance_cell)) |>
    count(concordance_cell, name = "claim_count") |>
    mutate(claim_proportion = claim_count / total_claims)
  
  study_ids <- sort(unique(claim_level$study_id))
  
  complete_grid <- expand.grid(
    study_id = study_ids, concordance_cell = concordance_cell_levels, stringsAsFactors = FALSE
  ) |> as_tibble()
  
  study_cell_counts <- claim_level |>
    mutate(concordance_cell = as.character(concordance_cell)) |>
    count(study_id, concordance_cell, name = "claims_in_cell")
  
  study_totals <- claim_level |> count(study_id, name = "claims_in_study")
  
  study_weighted <- complete_grid |>
    left_join(study_cell_counts, by = c("study_id", "concordance_cell")) |>
    mutate(claims_in_cell = coalesce(claims_in_cell, 0L)) |>
    left_join(study_totals, by = "study_id") |>
    mutate(within_study_proportion = claims_in_cell / claims_in_study) |>
    group_by(concordance_cell) |>
    summarise(
      study_weighted_proportion = mean(within_study_proportion),
      studies_with_cell = sum(claims_in_cell > 0),
      .groups = "drop"
    )
  
  tibble(concordance_cell = concordance_cell_levels) |>
    left_join(claim_counts, by = "concordance_cell") |>
    left_join(study_weighted, by = "concordance_cell") |>
    mutate(
      claim_count = coalesce(claim_count, 0L),
      claim_proportion = coalesce(claim_proportion, 0),
      study_weighted_proportion = coalesce(study_weighted_proportion, 0),
      studies_with_cell = coalesce(studies_with_cell, 0L),
      concordance_cell = factor(concordance_cell, levels = concordance_cell_levels)
    ) |>
    arrange(concordance_cell)
}

plot_evidence_plane <- function(claim_level, alpha, k, output_path) {
  p_threshold <- -log10(alpha)
  bf_threshold <- log10(k)
  
  plot <- ggplot(claim_level, aes(x = negative_log10_p, y = log10_bf10, colour = concordance_cell)) +
    geom_hline(yintercept = c(-bf_threshold, bf_threshold), linetype = "dashed") +
    geom_vline(xintercept = p_threshold, linetype = "dashed") +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_point(size = 2.7, alpha = 0.85) +
    labs(
      x = expression(-log[10](p)), y = expression(log[10](BF[10])), colour = "Concordance cell",
      title = "Frequentist-Bayesian evidence plane",
      subtitle = paste0(
        "Primary prior; alpha = ", alpha, "; H1 threshold BF10 >= ", k,
        "; H0 threshold BF10 <= ", round(1 / k, 3)
      )
    ) +
    theme_minimal(base_size = 12)
  
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = output_path, plot = plot, width = 8, height = 5.5, dpi = 300)
  invisible(plot)
}

build_detailed_rank_table <- function(claim_level) {
  claim_level |>
    select(
      claim_id, study_id, stat_test, p_value, p_band, frequentist_result, bf10, bf_strength,
      favoured_side, concordance_cell, concordance_status, negative_log10_p, log10_bf10,
      prior_sensitivity_span
    ) |>
    arrange(desc(log10_bf10))
}

plot_detailed_rank <- function(detailed_rank, output_path) {
  tr <- function(x) asinh(x)
  brks <- c(-5, -2, -1, 0, 1, 2, 5, 10, 30, 90, 220)
  df <- detailed_rank |> mutate(claim_id = fct_reorder(claim_id, log10_bf10))

  p <- ggplot(df, aes(x = tr(log10_bf10), y = claim_id)) +
    geom_vline(xintercept = tr(log10(c(1 / 3, 3))), linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
    geom_point(aes(colour = concordance_status), size = 3) +
    geom_text(aes(label = p_band), hjust = -0.15, size = 2.6, colour = "grey30") +
    scale_x_continuous(breaks = tr(brks), labels = brks) +
    scale_colour_manual(values = c("Concordant" = "#2b8a3e", "Inconclusive" = "#868e96", "Discordant" = "#e03131")) +
    labs(
      x = expression(log[10](BF[10]) ~ "(asinh scale; primary prior)"), y = NULL, colour = "Concordance status",
      title = "Detailed evidence rank",
      subtitle = paste0(
        "Ranked by Bayes-factor strength (Jeffreys 1961 categories); asinh axis keeps extreme and ",
        "near-threshold claims legible; label shows the p-value band; colour shows six-cell status"
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank())

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = output_path, plot = p, width = 9, height = max(4, 0.35 * nrow(df) + 1.5), dpi = 300)
  invisible(p)
}

plot_detailed_rank_zoom <- function(detailed_rank, output_path, lim = 0.8) {
  df <- detailed_rank |>
    filter(is.finite(log10_bf10), abs(log10_bf10) <= lim) |>
    mutate(claim_id = fct_reorder(claim_id, log10_bf10))

  p <- ggplot(df, aes(x = log10_bf10, y = claim_id)) +
    annotate("rect", xmin = log10(1 / 3), xmax = log10(3), ymin = -Inf, ymax = Inf, fill = "#9E9E9E", alpha = 0.12) +
    geom_vline(xintercept = c(log10(1 / 3), log10(3)), linetype = "dashed", linewidth = 0.35) +
    geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
    geom_point(aes(colour = concordance_status), size = 3) +
    geom_text(aes(label = p_band), hjust = -0.2, size = 2.6, colour = "grey30") +
    scale_colour_manual(values = c("Concordant" = "#2b8a3e", "Inconclusive" = "#868e96", "Discordant" = "#e03131")) +
    coord_cartesian(xlim = c(-lim - 0.15, lim + 0.15)) +
    labs(
      x = expression(log[10](BF[10]) ~ "(primary prior)"), y = NULL, colour = "Concordance status",
      title = "Detailed evidence rank — problematic zone",
      subtitle = paste0("Only claims with |log10 BF10| <= ", lim, ", where the six-cell decision is actually in play")
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank())

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = output_path, plot = p, width = 8, height = max(3, 0.42 * nrow(df) + 1.2), dpi = 300)
  invisible(p)
}

plot_problematic_pvalues <- function(detailed_rank, output_path, alpha = 0.05) {
  tr <- function(x) asinh(x)
  brks <- c(-2, -1, 0, 1, 2, 5, 10, 30, 90)
  fmt_p <- function(x) ifelse(x < 0.001, formatC(x, format = "e", digits = 1), formatC(x, format = "f", digits = 3))
  df <- detailed_rank |>
    filter(concordance_status %in% c("Inconclusive", "Discordant")) |>
    mutate(
      claim_id = fct_reorder(claim_id, log10_bf10),
      p_label = paste0("p = ", fmt_p(p_value), ifelse(p_value < alpha, " *", ""))
    )

  xr <- range(tr(df$log10_bf10))
  keep <- tr(brks) >= xr[1] - 0.3 & tr(brks) <= xr[2] + 0.3

  p <- ggplot(df, aes(x = tr(log10_bf10), y = claim_id)) +
    annotate("rect", xmin = tr(log10(1 / 3)), xmax = tr(log10(3)), ymin = -Inf, ymax = Inf, fill = "#9E9E9E", alpha = 0.12) +
    geom_vline(xintercept = tr(log10(c(1 / 3, 3))), linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
    geom_point(aes(colour = concordance_status), size = 3) +
    geom_text(aes(label = p_label), hjust = -0.15, size = 2.7, colour = "grey25") +
    scale_x_continuous(breaks = tr(brks)[keep], labels = brks[keep], expand = expansion(mult = c(0.05, 0.2))) +
    scale_colour_manual(values = c("Inconclusive" = "#868e96", "Discordant" = "#e03131")) +
    labs(
      x = expression(log[10](BF[10]) ~ "(asinh scale; primary prior)"), y = NULL, colour = "Concordance status",
      title = "Problematic claims: evidence and p-values",
      subtitle = paste0(
        "Only inconclusive and discordant claims; grey band = Bayesian inconclusive region; ",
        "label shows the frequentist p-value (* = significant at alpha = ", alpha, ")"
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank())

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = output_path, plot = p, width = 9, height = max(3, 0.4 * nrow(df) + 1.5), dpi = 300)
  invisible(p)
}

build_concordance_outputs <- function(
    results_path = "outputs/tables/bayes_factor_results.csv",
    alpha = .05,
    k = 3,
    claim_output_path = "outputs/tables/concordance_claim_level.csv",
    summary_output_path = "outputs/tables/concordance_summary.csv",
    figure_output_path = "outputs/figures/evidence_concordance_plane.png",
    detailed_rank_figure_path = "outputs/figures/detailed_evidence_rank.png",
    detailed_rank_zoom_figure_path = "outputs/figures/detailed_evidence_rank_zoom.png",
    problematic_pvalues_figure_path = "outputs/figures/problematic_claims_pvalues.png",
    remove_legacy_outputs = TRUE
) {
  if (remove_legacy_outputs) {
    legacy_paths <- c(
      "outputs/tables/agreement_rank.csv", "outputs/figures/agreement_rank.png",
      "outputs/figures/evidence_plane.png", "outputs/tables/prior_sensitivity.csv",
      "outputs/tables/detailed_evidence_rank.csv",
      "outputs/figures/evidence_concordance_plane.png"
    )
    unlink(legacy_paths[file.exists(legacy_paths)])
  }
  
  results <- read_csv(results_path, show_col_types = FALSE)
  classified <- add_evidence_classifications(results = results, alpha = alpha, k = k)
  prior_sensitivity <- summarise_prior_sensitivity(classified)
  
  claim_level <- classified |>
    filter(prior_label == "primary") |>
    left_join(prior_sensitivity, by = c("claim_id", "study_id")) |>
    select(
      claim_id, study_id, stat_test, p_value, p_band, frequentist_result, bf10, bf_conclusion,
      bf_strength, favoured_side, concordance_cell, concordance_status, negative_log10_p, log10_bf10,
      n_prior_specifications, minimum_bf10, maximum_bf10, prior_sensitivity_span,
      bf_conclusion_changed, favoured_side_changed
    ) |>
    arrange(study_id, claim_id)
  
  concordance_summary <- build_concordance_summary(claim_level)
  detailed_rank <- build_detailed_rank_table(claim_level)
  
  dir.create(dirname(claim_output_path), recursive = TRUE, showWarnings = FALSE)
  write_csv(claim_level, claim_output_path, na = "")
  write_csv(concordance_summary, summary_output_path, na = "")
  
  plot_detailed_rank(detailed_rank = detailed_rank, output_path = detailed_rank_figure_path)
  plot_detailed_rank_zoom(detailed_rank = detailed_rank, output_path = detailed_rank_zoom_figure_path)
  plot_problematic_pvalues(detailed_rank = detailed_rank, output_path = problematic_pvalues_figure_path, alpha = alpha)
  
  message(
    "Created concordance outputs for ", nrow(claim_level), " claims from ",
    n_distinct(claim_level$study_id), " studies."
  )
  
  invisible(list(
    claim_level = claim_level, summary = concordance_summary,
    prior_sensitivity = prior_sensitivity, detailed_rank = detailed_rank
  ))
}

if (sys.nframe() == 0L) {
  build_concordance_outputs()
}