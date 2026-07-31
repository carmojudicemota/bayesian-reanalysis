library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(ggrepel)
library(patchwork)

build_detailed_rank_table <- function(claim_level) {
  claim_level |>
    select(
      claim_id, study_id, stat_test, p_value, p_band, frequentist_result, bf10, bf_strength,
      favoured_side, concordance_cell, concordance_status, negative_log10_p, log10_bf10,
      prior_sensitivity_span
    ) |>
    arrange(desc(log10_bf10))
}
plot_concordace_squares <- function(
    in_path = "outputs/tables/concordance_summary.csv",
    out_path = "outputs/figures/concordance_squares.png") {
  library(dplyr); library(tidyr)
  d <- readr::read_csv(in_path, show_col_types = FALSE) |>
    tidyr::separate(concordance_cell, c("sig", "bf"), sep = " \\+ ", remove = FALSE) |>
    mutate(bf  = recode(bf, "inconclusive" = "Inconclusive"),
           sig = factor(sig, c("Nonsignificant", "Significant")),
           bf  = factor(bf,  c("H0", "Inconclusive", "H1")),
           status = case_when(
             concordance_cell %in% c("Significant + H1", "Nonsignificant + H0") ~ "Concordant",
             grepl("inconclusive", concordance_cell)                            ~ "Inconclusive",
             TRUE ~ "Discordant"))
  p <- ggplot(d, aes(bf, sig, fill = status, alpha = claim_proportion)) +
    geom_tile(colour = "white", linewidth = 3) +
    geom_text(aes(label = ifelse(claim_count > 0, paste0(claim_count, "\n", round(claim_proportion*100), "%"),"")),
              fontface = "bold", size = 5, alpha = 1)+
    scale_fill_manual(values = concordance_colours) +
    labs(x = "Bayes Factor Conclusion", y = NULL , fill = NULL,
         title = "3x2 Concordance Scheme") +
    guides(alpha = "none") +
    theme_reanalysis() + theme(panel.grid = element_blank(),
                               plot.title = element_text(hjust = 0.5))
  save_fig(p,out_path,w = 8, h = 4.5)
}
plot_evidence_plan <- function(
    in_path = "outputs/tables/concordance_claim_level.csv",
    out_path = "outputs/figures/evidence_plan.png",
    alpha = alpha_default, k = k_default) {
  library(dplyr)
  cl <- readr::read_csv(in_path, show_col_types = FALSE)
  sbb <- function(x) log10(-1 / (exp(1) * 10^(-x) * log(10^(-x))))
  p <- ggplot(cl, aes(negative_log10_p, log10_bf10, colour = concordance_status)) +
    annotate("rect", xmin=-Inf, xmax =Inf, ymin=log10(k), ymax=Inf,fill=tint_blue)+
    annotate("rect", xmin=-Inf, xmax =Inf, ymin=-Inf, ymax=-log10(k),fill=tint_red)+
    stat_function(fun = sbb, colour = blue_med, linetype = "dashed",
                  linewidth = 0.8, xlim= c(0.44,max(cl$negative_log10_p))) +
    geom_hline(yintercept = c(-log10(k),log10(k)), linetype = "dashed", colour = "grey40")+
    geom_vline(xintercept = -log10(alpha), linetype = "dashed", colour = "grey40")+
    geom_hline(yintercept = 0, linewidth = 0.3)+
    geom_point(size = 3) +
    scale_colour_manual(values = concordance_colours) +
    coord_cartesian(ylim = c(-1.3, 7.3)) +
    labs(x = expression(-log[10](p)), y = expression(log[10](BF[10])),
         colour = NULL, title = "Frequentist-Bayesian Evidence Plane") +
    theme_reanalysis()
  save_fig(p,out_path,w=8,h=5.5)
}
plot_detailed_rank <- function(
    in_path = "outputs/tables/concordance_claim_level.csv",
    out_path = "outputs/figures/detailed_rank.png") {
  library(dplyr); library(forcats)
  d <- readr::read_csv(in_path, show_col_types = FALSE) |>
    mutate(claim_id = fct_reorder(claim_id, log10_bf10))
  p <- ggplot(d, aes(log10_bf10, claim_id)) +
    geom_vline(xintercept = log10(jeffreys_bf), colour = "grey88", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
    geom_point(aes(colour = concordance_status), size = 3) +
    geom_text(aes(label = p_band), hjust = -0.15, size = 2.6, colour = "grey30") +
    scale_colour_manual(values = concordance_colours) +
    labs(x = expression(log[10](BF[10])~"(primary prior)"), y = NULL,
         colour = NULL, title = "Detailed Evidence Rank") +
    theme_reanalysis() + theme(panel.grid = element_blank())
  save_fig(p, out_path, w = 9, h = max(4, 0.35 * nrow(d) + 1.5))
}
plot_prior_sensitivity <-function(
    in_path = "outputs/tables/bayes_factor_results.csv",
    out_path = "outputs/figures/prior_sensitivity.png") {
  library(dplyr); library(tidyr); library(forcats)
  d <- readr::read_csv(in_path,show_col_types = FALSE) |>
    mutate(log10_bf10 = log10(bf10)) |>
    group_by(claim_id) |>
    mutate(primary_lbf = log10_bf10[prior_label == "primary"]) |>
    ungroup() |>
    mutate(claim_id = fct_reorder(claim_id, primary_lbf))
  spans <- d |> group_by(claim_id) |>
    summarise(lo = min(log10_bf10), hi = max(log10_bf10), .groups = "drop")
  p <- ggplot() +
    geom_vline(xintercept = c(log10(3), log10(1/3)), linetype = "dashed", colour = red_main) +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    geom_segment(data = spans, aes(x = lo, xend = hi, y = claim_id, yend = claim_id),
                 colour = "grey40", linewidth = 1.6) +
    geom_point(data = d, aes(log10_bf10, claim_id, colour = prior_label), size = 2.4) +
    scale_colour_manual(values = c(narrow = blue_light, primary = blue_deep, wide = red_main)) +
    labs(x = expression(log[10](BF[10])), y = NULL, colour = "Prior",
         title = "Prior Sensitivity (Narrow / Primary / Wide)") +
    coord_cartesian(xlim = c(-1.2, 8)) +
    theme_reanalysis()
  save_fig(p, out_path, w = 8, h = max(4, 0.35 * n_distinct(d$claim_id) + 1.5))
}
#inspired by the wetzels article
plot_evidence_grid <- function(
    in_path = "outputs/tables/concordance_claim_level.csv",
    out_path = "outputs/figures/evidence_grid_p.png") {
  library(dplyr)
  cl <- readr::read_csv(in_path, show_col_types = FALSE)
  p <- ggplot(cl, aes(p_value, bf10, colour = concordance_status)) +
    geom_hline(yintercept = jeffreys_bf, colour = "grey90", linewidth = 0.3) +
    geom_vline(xintercept = c(.001, .01, .05), colour = "grey90", linewidth = 0.3) +
    geom_point(size = 2.6) +
    scale_x_log10(breaks = c(1e-5, 1e-3, 1e-2, .05, 1),
                  limits = c(1e-6, 1), oob = scales::squish) +
    scale_y_log10(breaks = c(1/10, 1/3, 1, 3, 10, 30, 100),
                  labels = c("1/10","1/3","1","3","10","30","100"),
                  limits = c(0.1, 1e5), oob = scales::squish) + scale_colour_manual(values = concordance_colours) +
    labs(x = "p-value", y = "Bayes Factor", colour = NULL,
         title = "Evidence Grid: Bayes Factor vs P",
         subtitle = "Jeffreys categories; after Wetzels et al. (2011)") +
    theme_reanalysis()
  save_fig(p, out_path, w = 8, h = 5.5)
}
add_margin_pct <- function(p, cl) {
  bands <- cut(cl$p_value, c(-1,.001,.01,.05,1),
               labels = c("<.001",".001-.01",".01-.05",">.05"))
  tab <- round(100 * prop.table(table(bands)))
  p + annotate("text", x = c(3e-4,3e-3,.022,.35), y = 1.4e5,
               label = paste0(tab, "%"), size = 3, colour = "grey35")
  
}





status_cols <- c(Concordant = "#2166AC", Inconclusive = "#9E9E9E", Discordant = "#B2182B")
prior_cols  <- c(narrow = "#92C5DE", primary = "#08306B", wide = "#B2182B")
thr_p <- -log10(0.05)
thr_b <- log10(3)
e_const <- exp(1)

theme_project <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

short_claim <- function(x) gsub("_claim_", "-c", gsub("study_", "S", x))
safe_numeric <- function(x) suppressWarnings(as.numeric(x))

sbb_bf01_min <- function(p) ifelse(p > 0 & p < 1 / e_const, -e_const * p * log(p), NA_real_)
sbb_max_log10_bf10 <- function(p) { b <- sbb_bf01_min(p); ifelse(is.na(b), NA_real_, -log10(b)) }
sbb_min_post_h0 <- function(p, prior_odds_h0 = 1) {
  b <- sbb_bf01_min(p) * prior_odds_h0
  ifelse(is.na(b), NA_real_, b / (1 + b))
}

cl <- read_csv("outputs/tables/concordance_claim_level.csv", show_col_types = FALSE)
bf <- read_csv("outputs/tables/bayes_factor_results.csv", show_col_types = FALSE)

if (!"family" %in% names(cl)) {
  family_lookup <- bf |>
    transmute(
      claim_id,
      family = coalesce(
        if ("stat_test" %in% names(bf)) as.character(stat_test) else NA_character_,
        if ("bf_family" %in% names(bf)) as.character(bf_family) else NA_character_,
        if ("method" %in% names(bf)) as.character(method) else NA_character_
      )
    ) |>
    filter(!is.na(family)) |>
    distinct(claim_id, .keep_all = TRUE)
  cl <- cl |> left_join(family_lookup, by = "claim_id")
}

cl <- cl |>
  mutate(
    negative_log10_p = safe_numeric(negative_log10_p),
    log10_bf10 = safe_numeric(log10_bf10),
    p_value = safe_numeric(p_value),
    concordance_status = factor(
      concordance_status,
      levels = c("Concordant", "Inconclusive", "Discordant")
    )
  )

bf <- bf |>
  mutate(log10_bf10 = safe_numeric(log10_bf10), prior_label = as.character(prior_label))

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

zoom_lim <- 1.5
cat_of <- function(x) ifelse(x > thr_b, "H1", ifelse(x < -thr_b, "H0", "Inconclusive"))
tr_asinh <- function(x) asinh(x)

draw_dumbbell <- function(df, xlim, title, subtitle) {
  df <- df |> dplyr::mutate(row = rank(primary, ties.method = "first"))
  ggplot(df) +
    annotate("rect", xmin = -thr_b, xmax = thr_b, ymin = -Inf, ymax = Inf, fill = "#9E9E9E", alpha = 0.12) +
    geom_vline(xintercept = c(-thr_b, thr_b), linetype = "dashed", linewidth = 0.35) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_segment(aes(x = lo, xend = hi, y = row, yend = row, colour = connector), linewidth = 1.4) +
    geom_point(aes(nar, row), colour = prior_cols[["narrow"]],  size = 2.1, na.rm = TRUE) +
    geom_point(aes(pri, row), colour = prior_cols[["primary"]], size = 3.0, na.rm = TRUE) +
    geom_point(aes(wid, row), colour = prior_cols[["wide"]],    size = 2.1, na.rm = TRUE) +
    scale_colour_manual(values = c(changes = "#D55E00", stable = "#4477AA"), name = NULL,
                        labels = c(changes = "Category changes", stable = "Category stable")) +
    scale_y_continuous(breaks = df$row, labels = short_claim(df$claim_id), expand = expansion(add = 0.8)) +
    coord_cartesian(xlim = xlim) +
    labs(title = title, subtitle = subtitle, x = expression(log[10](BF[10])), y = NULL) +
    theme_project() + theme(axis.text.y = element_text(size = 7))
}

# ---- Evidence plane -----------------------------------------
lab_df <- cl |> filter(concordance_status == "Inconclusive")

evidence_plane_main <- ggplot(cl, aes(negative_log10_p, log10_bf10)) +
  annotate("rect", xmin = thr_p, xmax = Inf, ymin = thr_b, ymax = Inf, fill = "#009E73", alpha = 0.06) +
  annotate("rect", xmin = -Inf, xmax = thr_p, ymin = -Inf, ymax = -thr_b, fill = "#009E73", alpha = 0.06) +
  annotate("rect", xmin = thr_p, xmax = Inf, ymin = -Inf, ymax = -thr_b, fill = "#D55E00", alpha = 0.07) +
  annotate("rect", xmin = -Inf, xmax = thr_p, ymin = thr_b, ymax = Inf, fill = "#D55E00", alpha = 0.07) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -thr_b, ymax = thr_b, fill = "#999999", alpha = 0.10) +
  geom_hline(yintercept = c(-thr_b, thr_b), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = thr_p, linetype = "dashed", linewidth = 0.4) +
  geom_point(aes(colour = concordance_status), size = 2.8) +
  geom_text_repel(
    data = lab_df,
    aes(label = short_claim(claim_id)),
    size = 2.8, min.segment.length = 0, max.overlaps = Inf
  ) +
  scale_colour_manual(values = status_cols, name = NULL, drop = FALSE) +
  coord_cartesian(xlim = c(0, 5.2), ylim = c(-3.2, 3.2)) +
  labs(
    title = "Frequentist–Bayesian Evidence Plane",
    subtitle = "Decision zone; extreme Bayes factors are clipped at the panel edge",
    x = expression(Frequentist~evidence~~-log[10](p)),
    y = expression(Bayesian~evidence~~log[10](BF[10]))
  ) +
  theme_project()

evidence_plane_strip <- ggplot(cl, aes(x = 0, y = log10_bf10, colour = concordance_status)) +
  annotate("rect", xmin = -0.6, xmax = 0.6, ymin = -thr_b, ymax = thr_b, fill = "#999999", alpha = 0.12) +
  geom_jitter(width = 0.28, height = 0, size = 1.8, show.legend = FALSE) +
  scale_colour_manual(values = status_cols, drop = FALSE) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(x = NULL, y = expression(log[10](BF[10])~"(symlog)"), title = "Full spread") +
  theme_project() +
  theme(axis.text.x = element_blank())

evidence_plane <- evidence_plane_main + evidence_plane_strip + plot_layout(widths = c(3, 1))
ggsave("outputs/figures/evidence_plane.png", evidence_plane, width = 11, height = 7, dpi = 300)


grid <- expand_grid(freq = c("Significant", "Nonsignificant"),bayes = c("H0", "Inconclusive", "H1")) |>
  mutate(
    status = case_when(
      freq == "Significant" & bayes == "H1" ~ "Concordant",
      freq == "Nonsignificant" & bayes == "H0" ~ "Concordant",
      bayes == "Inconclusive" ~ "Inconclusive",
      TRUE ~ "Discordant"
    ),
    x = match(bayes, c("H0", "Inconclusive", "H1")),
    y = match(freq, c("Nonsignificant", "Significant"))
  )

set.seed(123)
pts <- cl |>
  mutate(
    freq = if_else(frequentist_result == "Significant", "Significant", "Nonsignificant"),
    bayes = as.character(bf_conclusion),
    x = match(bayes, c("H0", "Inconclusive", "H1")),
    y = match(freq, c("Nonsignificant", "Significant")),
    xj = x + runif(n(), -0.27, 0.27),
    yj = y + runif(n(), -0.23, 0.08)
  ) |>
  left_join(grid |> select(freq, bayes, status), by = c("freq", "bayes"))

counts <- cl |>
  mutate(
    freq = if_else(frequentist_result == "Significant", "Significant", "Nonsignificant"),
    bayes = as.character(bf_conclusion)
  ) |>
  count(freq, bayes, name = "n") |>
  mutate(
    x = match(bayes, c("H0", "Inconclusive", "H1")),
    y = match(freq, c("Nonsignificant", "Significant"))
  )

concordance_plot <- ggplot(grid, aes(x, y)) +
  geom_tile(aes(fill = status), alpha = 0.18, colour = "grey85", width = 0.98, height = 0.98) +
  geom_text(data = counts, aes(x = x - 0.32, y = y + 0.30, label = n),
            inherit.aes = FALSE, size = 7, fontface = "bold") +
  geom_point(data = pts, aes(x = xj, y = yj, colour = status),
             inherit.aes = FALSE, size = 2.4) +
  scale_fill_manual(values = status_cols, guide = "none", drop = FALSE) +
  scale_colour_manual(values = status_cols, guide = "none", drop = FALSE) +
  scale_x_continuous(breaks = 1:3,
                     labels = c("Bayesian: H0", "Bayesian: Inconclusive", "Bayesian: H1"),
                     position = "top", expand = expansion(add = 0.5)) +
  scale_y_continuous(breaks = 1:2, labels = c("Nonsignificant", "Significant"),
                     expand = expansion(add = 0.5)) +
  coord_fixed() +
  labs(title = "Six-Cell Concordance Matrix",
       subtitle = "Frequentist vs Bayesian Evidence", x = NULL, y = NULL) +
  theme_project()

ggsave("outputs/figures/concordance.png", concordance_plot, width = 9, height = 6, dpi = 300)


bars_df <- cl |>
  filter(is.finite(log10_bf10)) |>
  mutate(claim = reorder(short_claim(claim_id), log10_bf10))

rank_brks <- c(-5, -2, -1, 0, 1, 2, 5, 10, 30, 90, 220)

evidence_bars <- ggplot(bars_df, aes(claim, tr_asinh(log10_bf10), fill = log10_bf10)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = tr_asinh(-thr_b), ymax = tr_asinh(thr_b),
           fill = "#9E9E9E", alpha = 0.12) +
  geom_hline(yintercept = c(tr_asinh(-thr_b), tr_asinh(thr_b)), linetype = "dashed", linewidth = 0.35) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_y_continuous(breaks = tr_asinh(rank_brks), labels = rank_brks) +
  scale_fill_gradient2(low = "#B2182B", mid = "#EAEAEA", high = "#2166AC", midpoint = 0,
                       name = expression(log[10](BF[10]))) +
  labs(title = "Detailed Evidence Rank",
       subtitle = "asinh axis so extreme and near-threshold claims are both legible", x = NULL,
       y = expression(log[10](BF[10])~"(asinh)")) +
  theme_project() +
  theme(axis.text.y = element_text(size = 7))

ggsave("outputs/figures/diverging_evidence_bars.png", evidence_bars,
       width = 9, height = 4 + 0.20 * nrow(bars_df), dpi = 300, limitsize = FALSE)

bars_zoom_df <- cl |>
  filter(is.finite(log10_bf10), abs(log10_bf10) <= zoom_lim) |>
  mutate(claim = reorder(short_claim(claim_id), log10_bf10))

evidence_bars_zoom <- ggplot(bars_zoom_df, aes(claim, log10_bf10, fill = log10_bf10)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -thr_b, ymax = thr_b, fill = "#9E9E9E", alpha = 0.12) +
  geom_hline(yintercept = c(-thr_b, thr_b), linetype = "dashed", linewidth = 0.35) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_gradient2(low = "#B2182B", mid = "#EAEAEA", high = "#2166AC", midpoint = 0,
                       name = expression(log[10](BF[10]))) +
  labs(title = "Detailed Evidence Rank — threshold zoom",
       subtitle = paste0("Claims with |log10 BF10| ≤ ", zoom_lim), x = NULL,
       y = expression(log[10](BF[10]))) +
  theme_project() +
  theme(axis.text.y = element_text(size = 8))

ggsave("outputs/figures/diverging_evidence_bars_zoom.png", evidence_bars_zoom,
       width = 8, height = 3 + 0.30 * nrow(bars_zoom_df), dpi = 300, limitsize = FALSE)

prior_variants <- bf |>
  dplyr::filter(is.finite(log10_bf10)) |>
  dplyr::distinct(claim_id, prior_label, log10_bf10)

prior_summ <- prior_variants |>
  dplyr::group_by(claim_id) |>
  dplyr::summarise(
    primary = log10_bf10[match("primary", prior_label)],
    lower = min(log10_bf10), upper = max(log10_bf10),
    category_changes = dplyr::n_distinct(cat_of(log10_bf10)) > 1,
    .groups = "drop"
  ) |>
  dplyr::mutate(primary = ifelse(is.na(primary), (lower + upper) / 2, primary))

prior_nw <- prior_variants |>
  dplyr::filter(prior_label %in% c("narrow", "wide")) |>
  tidyr::pivot_wider(names_from = prior_label, values_from = log10_bf10)

prior_wide <- prior_summ |>
  dplyr::left_join(prior_nw, by = "claim_id") |>
  dplyr::mutate(connector = dplyr::if_else(category_changes, "changes", "stable"))
if (!"narrow" %in% names(prior_wide)) prior_wide$narrow <- NA_real_
if (!"wide"   %in% names(prior_wide)) prior_wide$wide   <- NA_real_

full_df <- prior_wide |>
  dplyr::transmute(claim_id, connector, primary,
                   lo = pmax(lower, -6), hi = pmin(upper, 6),
                   nar = pmax(pmin(narrow, 6), -6), pri = pmax(pmin(primary, 6), -6),
                   wid = pmax(pmin(wide, 6), -6))

prior_robustness <- draw_dumbbell(full_df, c(-6.4, 6.4),
  "Prior Robustness (all claims)",
  "Single point = primary-only claim; axis capped at ±6")
ggsave("outputs/figures/prior_robustness.png", prior_robustness,
       width = 9, height = 4 + 0.22 * nrow(full_df), dpi = 300, limitsize = FALSE)

zoom_df <- prior_wide |>
  dplyr::filter(abs(primary) <= zoom_lim) |>
  dplyr::transmute(claim_id, connector, primary,
                   lo = lower, hi = upper, nar = narrow, pri = primary, wid = wide)

prior_robustness_zoom <- draw_dumbbell(zoom_df, c(-zoom_lim - 0.4, zoom_lim + 0.4),
  "Prior Robustness — threshold zoom",
  paste0("Claims with |log10 BF10| ≤ ", zoom_lim, "; full span, no capping"))
ggsave("outputs/figures/prior_robustness_zoom.png", prior_robustness_zoom,
       width = 9, height = 4 + 0.30 * nrow(zoom_df), dpi = 300, limitsize = FALSE)


prior_long <- bf |>
  filter(prior_label %in% c("narrow", "primary", "wide")) |>
  select(claim_id, prior_label, log10_bf10) |>
  distinct() |>
  group_by(claim_id) |>
  filter(n_distinct(prior_label) >= 2) |>
  ungroup() |>
  mutate(
    prior_label = factor(prior_label, levels = c("narrow", "primary", "wide")),
    capped = pmax(pmin(log10_bf10, 6), -6),
    category = case_when(log10_bf10 > thr_b ~ "H1", log10_bf10 < -thr_b ~ "H0", TRUE ~ "Inconclusive")
  )

prior_long <- prior_long |>
  left_join(
    prior_long |> group_by(claim_id) |> summarise(category_changes = n_distinct(category) > 1, .groups = "drop"),
    by = "claim_id"
  )

prior_slope_graph <- ggplot(prior_long, aes(prior_label, capped, group = claim_id)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -thr_b, ymax = thr_b, fill = "#9E9E9E", alpha = 0.12) +
  geom_hline(yintercept = c(-thr_b, thr_b), linetype = "dashed", linewidth = 0.35) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(aes(colour = category_changes), linewidth = 0.8, alpha = 0.75) +
  geom_point(aes(fill = prior_label), shape = 21, colour = "white", stroke = 0.4, size = 2.5) +
  scale_colour_manual(values = c(`FALSE` = "#4477AA", `TRUE` = "#D55E00"), name = NULL,
                      labels = c(`FALSE` = "Category stable", `TRUE` = "Category changes")) +
  scale_fill_manual(values = prior_cols, name = "Prior") +
  scale_x_discrete(labels = c(narrow = "Narrow", primary = "Primary", wide = "Wide")) +
  coord_cartesian(ylim = c(-6.4, 6.4)) +
  labs(title = "Prior Slope Graph",
       subtitle = "Movement Across the Prior Grid; Axis Capped at ±6",
       x = NULL, y = expression(log[10](BF[10]))) +
  theme_project()

ggsave("outputs/figures/prior_slope_graph.png", prior_slope_graph, width = 9, height = 7, dpi = 300)

prim_lookup <- prior_long |> filter(prior_label == "primary") |> distinct(claim_id, prim = log10_bf10)
slope_zoom_df <- prior_long |> left_join(prim_lookup, by = "claim_id") |> filter(abs(prim) <= zoom_lim)

prior_slope_graph_zoom <- ggplot(slope_zoom_df, aes(prior_label, log10_bf10, group = claim_id)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -thr_b, ymax = thr_b, fill = "#9E9E9E", alpha = 0.12) +
  geom_hline(yintercept = c(-thr_b, thr_b), linetype = "dashed", linewidth = 0.35) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(aes(colour = category_changes), linewidth = 0.8, alpha = 0.75) +
  geom_point(aes(fill = prior_label), shape = 21, colour = "white", stroke = 0.4, size = 2.5) +
  scale_colour_manual(values = c(`FALSE` = "#4477AA", `TRUE` = "#D55E00"), name = NULL,
                      labels = c(`FALSE` = "Category stable", `TRUE` = "Category changes")) +
  scale_fill_manual(values = prior_cols, name = "Prior") +
  scale_x_discrete(labels = c(narrow = "Narrow", primary = "Primary", wide = "Wide")) +
  coord_cartesian(ylim = c(-zoom_lim - 0.4, zoom_lim + 0.4)) +
  labs(title = "Prior Slope Graph — threshold zoom",
       subtitle = paste0("Claims with |log10 BF10| ≤ ", zoom_lim, "; no capping"),
       x = NULL, y = expression(log[10](BF[10]))) +
  theme_project()

ggsave("outputs/figures/prior_slope_graph_zoom.png", prior_slope_graph_zoom, width = 9, height = 6, dpi = 300)


grades <- c("\u22641/100", "1/100–1/10", "1/10–1/3", "1/3–3", "3–10", "10–100", "\u2265100")
grade_cols <- c(
  "\u22641/100" = "#67001F", "1/100–1/10" = "#B2182B", "1/10–1/3" = "#F4A582",
  "1/3–3" = "#EAEAEA", "3–10" = "#92C5DE", "10–100" = "#4393C3", "\u2265100" = "#2166AC"
)
grade_df <- cl |>
  filter(!is.na(family), is.finite(log10_bf10)) |>
  mutate(grade = cut(log10_bf10, breaks = c(-Inf, -2, -1, -thr_b, thr_b, 1, 2, Inf),
                     labels = grades, include.lowest = TRUE, ordered_result = TRUE)) |>
  count(family, grade, name = "n", .drop = FALSE) |>
  complete(family, grade = factor(grades, levels = grades, ordered = TRUE), fill = list(n = 0)) |>
  group_by(family) |>
  mutate(
    total = sum(n),
    share = if_else(total > 0, n / total, 0),
    negative = grade %in% grades[1:3],
    inconclusive = grade == "1/3–3",
    positive = grade %in% grades[5:7],
    neg_total = sum(share[negative]),
    inc_share = sum(share[inconclusive]),
    start = -neg_total - inc_share / 2,
    order_index = as.integer(grade)
  ) |>
  arrange(family, order_index) |>
  mutate(xmin = start + lag(cumsum(share), default = 0), xmax = start + cumsum(share)) |>
  ungroup()

family_order <- grade_df |>
  group_by(family) |>
  summarise(positive_share = sum(share[positive]), negative_share = sum(share[negative]), .groups = "drop") |>
  arrange(positive_share - negative_share) |>
  pull(family)

grade_df <- grade_df |> mutate(family = factor(family, levels = family_order))

diverging_jeffreys_by_family <- ggplot(grade_df) +
  geom_rect(aes(xmin = xmin, xmax = xmax,
                ymin = as.numeric(family) - 0.38, ymax = as.numeric(family) + 0.38, fill = grade)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = grade_cols, drop = FALSE, name = "Jeffreys grade") +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1), limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  scale_y_continuous(breaks = seq_along(family_order), labels = family_order, expand = expansion(add = 0.6)) +
  labs(title = "Diverging Jeffreys Grades by Family",
       subtitle = "The inconclusive grade is centred on zero", x = "Share of Claims", y = NULL) +
  theme_project() +
  theme(axis.text.y = element_text(size = 8))

ggsave("outputs/figures/diverging_jeffreys_grades_by_family.png", diverging_jeffreys_by_family,
       width = 10, height = 4 + 0.35 * length(family_order), dpi = 300, limitsize = FALSE)

invisible(list(
  evidence_plane = evidence_plane,
  concordance = concordance_plot,
  diverging_evidence_bars = evidence_bars,
  diverging_evidence_bars_zoom = evidence_bars_zoom,
  prior_robustness = prior_robustness,
  prior_robustness_zoom = prior_robustness_zoom,
  prior_slope_graph = prior_slope_graph,
  prior_slope_graph_zoom = prior_slope_graph_zoom,
  diverging_jeffreys_grades_by_family = diverging_jeffreys_by_family
))


build_all_figures <- function() {
  plot_evidence_plan()
  plot_concordace_squares()
  plot_prior_sensitivity()
  plot_evidence_grid()
}