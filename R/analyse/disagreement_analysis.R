library(dplyr)
library(readr)
library(tidyr)

bf_path <- "outputs/tables/bayes_factor_results.csv"
claims_path <- "data/derived/claims.csv"
output_dir <- "outputs/tables"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bf <- read_csv(bf_path, show_col_types = FALSE)
claims <- read_csv(claims_path, show_col_types = FALSE)


classify_bf <- function(bf10) {
  case_when(bf10 >= 3   ~ "H1", bf10 <= 1/3 ~ "H0", TRUE~ "Inconclusive")
}

primary <- bf %>%
  filter(prior_label == "primary") %>%
  left_join(
    claims %>%
      select(
        claim_id,
        role,
        effect_size_type,
        effect_size_value,
        estimate,
        se_estimate,
        n_eff
      ),
    by = "claim_id"
  ) %>%
  mutate(
    frequentist_class = if_else(p_value < 0.05, "Significant","Non-significant"),
    bayesian_class = classify_bf(bf10),
    result_group = case_when(
      frequentist_class == "Significant" & bayesian_class == "H1" ~ "Significant + H1",
      frequentist_class == "Significant" & bayesian_class == "Inconclusive" ~ "Significant + inconclusive",
      frequentist_class == "Significant" & bayesian_class == "H0" ~ "Significant + H0",
      frequentist_class == "Non-significant" & bayesian_class == "H0" ~ "Non-significant + H0",
      frequentist_class == "Non-significant" & bayesian_class == "Inconclusive" ~ "Non-significant + inconclusive",
      frequentist_class == "Non-significant" & bayesian_class == "H1" ~ "Non-significant + H1")
  )

sensitivity <- bf %>%
  mutate(evidence_class = classify_bf(bf10), log10_bf10_calc = log10(bf10)) %>%
  group_by(claim_id) %>%
  summarise(
    n_specifications = n_distinct(prior_label),
    log10_bf_span = max(log10_bf10_calc, na.rm = TRUE) - min(log10_bf10_calc, na.rm = TRUE),
    n_evidence_classes = n_distinct(evidence_class),
    class_changed = n_evidence_classes > 1,
    .groups = "drop"
  )

analysis_data <- primary %>%
  left_join(sensitivity, by = "claim_id")

disagreement_drivers_summary <- analysis_data %>%
  group_by(result_group) %>%
  summarise(
    n_claims = n(),
    median_p = median(p_value, na.rm = TRUE),
    median_n = median(n_total, na.rm = TRUE),
    median_bf10 = median(bf10, na.rm = TRUE),
    sensitivity_evaluable = sum(n_specifications >= 3, na.rm = TRUE),
    class_changes = sum(class_changed & n_specifications >= 3, na.rm = TRUE),
    class_change_percent = if_else(sensitivity_evaluable > 0, 100 * class_changes / sensitivity_evaluable, NA_real_ ),
    median_log10_bf_span = if_else(
      sensitivity_evaluable > 0,
      median(log10_bf_span[n_specifications >= 3], na.rm = TRUE),
      NA_real_
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    result_group = factor(
      result_group,
      levels = c(
        "Significant + H1",
        "Significant + inconclusive",
        "Significant + H0",
        "Non-significant + H0",
        "Non-significant + inconclusive",
        "Non-significant + H1"
      )
    )
  ) %>%
  arrange(result_group) %>%
  filter(n_claims > 0) %>%
  mutate(
    median_p = signif(median_p, 4),
    median_n = round(median_n, 1),
    median_bf10 = signif(median_bf10, 4),
    class_change_percent = round(class_change_percent, 1),
    median_log10_bf_span = round(median_log10_bf_span, 3)
  )

write_csv(disagreement_drivers_summary, file.path(output_dir, "disagreement_drivers_summary.csv"))

print(disagreement_drivers_summary)


significant_p_band_summary <- primary %>%
  filter(frequentist_class == "Significant") %>%
  mutate(
    p_band = case_when(
      p_value < 0.005 ~ "p < .005",
      p_value < 0.010 ~ ".005 <= p < .01",
      p_value < 0.020 ~ ".01 <= p < .02",
      p_value < 0.030 ~ ".02 <= p < .03",
      TRUE ~ ".03 <= p < .05"
    ),
    p_band = factor(
      p_band,
      levels = c("p < .005",".005 <= p < .01",".01 <= p < .02",".02 <= p < .03",".03 <= p < .05"))
  ) %>%
  group_by(p_band) %>%
  summarise(
    n_claims = n(),
    bayesian_h1 = sum(bayesian_class == "H1"),
    bayesian_inconclusive =
      sum(bayesian_class == "Inconclusive"),
    bayesian_h0 = sum(bayesian_class == "H0"),
    inconclusive_percent =
      100 * bayesian_inconclusive / n_claims,
    .groups = "drop"
  ) %>%
  arrange(p_band) %>%
  mutate(inconclusive_percent = round(inconclusive_percent, 1))

write_csv(significant_p_band_summary, file.path(output_dir, "significant_p_band_summary.csv"))

print(significant_p_band_summary)

