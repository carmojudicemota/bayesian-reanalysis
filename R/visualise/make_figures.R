library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(ggrepel)
library(patchwork)

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
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(2, 2, 2, 2)
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
  mutate(bf10 = safe_numeric(bf10),
         log10_bf10 = dplyr::coalesce(safe_numeric(log10_bf10), log10(bf10)),
         prior_label = as.character(prior_label))

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

zoom_lim <- 1.5
cat_of <- function(x) ifelse(x > thr_b, "H1", ifelse(x < -thr_b, "H0", "Inconclusive"))
tr_asinh <- function(x) asinh(x)

jeff_bf10   <- c(1/100, 1/30, 1/10, 1/3, 1, 3, 10, 30, 100, 1e4, 1e10, 1e30, 1e90)
jeff_log10  <- log10(jeff_bf10)
jeff_labels <- c("1/100", "1/30", "1/10", "1/3", "1", "3", "10", "30", "100",
                 "10^4", "10^10", "10^30", "10^90")
jeff_lines  <- log10(c(1/100, 1/30, 1/10, 1/3, 3, 10, 30, 100))
jeff_bands  <- tibble::tibble(
  from = c(-7, -2, -1, -thr_b, thr_b, 1, 2),
  to   = c(-2, -1, -thr_b, thr_b, 1, 2, 230),
  band = factor(c("Decisive H0", "Strong H0", "Moderate H0", "Inconclusive",
                  "Moderate H1", "Strong H1", "Decisive H1"),
                levels = c("Decisive H0", "Strong H0", "Moderate H0", "Inconclusive",
                           "Moderate H1", "Strong H1", "Decisive H1"))
)
band_fill <- c("Decisive H0" = "#B2182B", "Strong H0" = "#E68A7F", "Moderate H0" = "#FADBD5",
               "Inconclusive" = "#ECECEC", "Moderate H1" = "#CFE0F1", "Strong H1" = "#7FB0D8",
               "Decisive H1" = "#2166AC")
jeff_sec_axis <- function() ggplot2::sec_axis(~ ., breaks = tr_asinh(jeff_log10),
                                              labels = sprintf("%.2f", jeff_log10),
                                              name = expression(log[10](BF[10])))

draw_dumbbell <- function(df, xlim, title, subtitle) {
  df <- df |> dplyr::mutate(row = rank(primary, ties.method = "first"))
  flag <- df |> dplyr::filter(claim_id == "study_40_claim_01")
  keep <- jeff_log10 >= xlim[1] & jeff_log10 <= xlim[2]
  vl <- jeff_lines[jeff_lines >= xlim[1] & jeff_lines <= xlim[2]]
  ggplot(df) +
    geom_rect(data = jeff_bands, aes(xmin = from, xmax = to, ymin = -Inf, ymax = Inf, fill = band),
              inherit.aes = FALSE, alpha = 0.45) +
    geom_vline(xintercept = vl, colour = "white", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
    geom_segment(aes(x = lo, xend = hi, y = row, yend = row, colour = connector), linewidth = 1.4) +
    geom_point(aes(nar, row), colour = prior_cols[["narrow"]],  size = 2.1, na.rm = TRUE) +
    geom_point(aes(pri, row), colour = prior_cols[["primary"]], size = 3.0, na.rm = TRUE) +
    geom_point(aes(wid, row), colour = prior_cols[["wide"]],    size = 2.1, na.rm = TRUE) +
    geom_segment(data = flag, aes(x = lo, xend = hi, y = row, yend = row), colour = "grey20", linewidth = 0.6, linetype = "22", inherit.aes = FALSE) +
    geom_point(data = flag, aes(pri, row), shape = 21, size = 3.8, stroke = 1.1, colour = "grey15", fill = NA, inherit.aes = FALSE) +
    geom_text(data = flag, aes(pri, row + 0.45, label = "classified by WAIC/LOO (BF prior-sensitive)"), size = 2.3, colour = "grey15", fontface = "italic", inherit.aes = FALSE) +
    scale_fill_manual(values = band_fill, name = "Jeffreys grade") +
    scale_colour_manual(values = c(changes = "#D55E00", stable = "#4477AA"), name = NULL,
                        labels = c(changes = "Category changes", stable = "Category stable")) +
    scale_x_continuous(breaks = jeff_log10[keep], labels = jeff_labels[keep],
                       sec.axis = sec_axis(~ ., breaks = jeff_log10[keep],
                                           labels = sprintf("%.2f", jeff_log10[keep]),
                                           name = expression(log[10](BF[10])))) +
    scale_y_continuous(breaks = df$row, labels = short_claim(df$claim_id), expand = expansion(add = 0.8)) +
    coord_cartesian(xlim = xlim) +
    guides(fill = guide_legend(nrow = 2, order = 1),
           colour = guide_legend(nrow = 1, order = 2)) +
    labs(title = title, subtitle = subtitle, x = expression(BF[10]~"(Jeffreys scale)"), y = NULL) +
    theme_project() + theme(axis.text.y = element_text(size = 7), panel.grid.major.x = element_blank())
}

lab_df <- cl |> filter(concordance_status == "Inconclusive")

evidence_plane_main <- ggplot(cl, aes(negative_log10_p, log10_bf10)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = thr_b, ymax = Inf, fill = "#2166AC", alpha = 0.06) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -thr_b, fill = "#B2182B", alpha = 0.06) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -thr_b, ymax = thr_b, fill = "#999999", alpha = 0.10) +
  geom_hline(yintercept = c(-thr_b, thr_b), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = c(-2, -1, 1, 2), linetype = "dotted", linewidth = 0.3, colour = "grey55") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = thr_p, linetype = "dashed", linewidth = 0.4) +
  geom_point(aes(colour = concordance_status), size = 2.8) +
  geom_text_repel(
    data = lab_df,
    aes(label = short_claim(claim_id)),
    size = 2.8, min.segment.length = 0, max.overlaps = Inf
  ) +
  scale_colour_manual(values = status_cols, name = NULL, drop = FALSE) +
  scale_y_continuous(
    breaks = c(-2, -1, -thr_b, 0, thr_b, 1, 2),
    labels = c("1/100", "1/10", "1/3", "1", "3", "10", "100"),
    sec.axis = sec_axis(~ ., breaks = c(-2, -1, -thr_b, 0, thr_b, 1, 2),
                        labels = sprintf("%.2f", c(-2, -1, -thr_b, 0, thr_b, 1, 2)),
                        name = expression(log[10](BF[10])))
  ) +
  scale_x_continuous(
    breaks = -log10(c(.05, .01, .001, .0001)),
    labels = sprintf("%.1f", -log10(c(.05, .01, .001, .0001))),
    sec.axis = sec_axis(~ ., breaks = -log10(c(.05, .01, .001, .0001)),
                        labels = c(".05", ".01", ".001", ".0001"), name = "p-value")
  ) +
  coord_cartesian(xlim = c(0, 5.2), ylim = c(-2.2, 2.2)) +
  labs(
    title = "Frequentist–Bayesian Evidence Plane",
    subtitle = "Left axis = BF10 on the Jeffreys scale; dashed = decision thresholds; zoomed to the ±100 decision zone",
    x = expression(Frequentist~evidence~~-log[10](p)),
    y = expression(BF[10])
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

counts_raw <- cl |>
  mutate(
    freq = if_else(frequentist_result == "Significant", "Significant", "Nonsignificant"),
    bayes = as.character(bf_conclusion)
  ) |>
  count(freq, bayes, name = "n")

counts <- grid |>
  transmute(freq, bayes, x, y, status) |>
  left_join(counts_raw, by = c("freq", "bayes")) |>
  mutate(
    n = coalesce(n, 0L),
    pct = 100 * n / sum(n),
    lab = paste0(n, "\n", sprintf("%.0f%%", pct))
  )

concordance_plot <- ggplot(grid, aes(x, y)) +
  geom_tile(aes(fill = status), alpha = 0.20, colour = "grey85", width = 0.98, height = 0.98) +
  geom_text(data = counts, aes(x = x, y = y + 0.29, label = lab),
            inherit.aes = FALSE, size = 5, fontface = "bold", lineheight = 0.9) +
  geom_point(data = pts, aes(x = xj, y = yj, colour = status),
             inherit.aes = FALSE, size = 2.2) +
  scale_fill_manual(values = status_cols, guide = "none", drop = FALSE) +
  scale_colour_manual(values = status_cols, guide = "none", drop = FALSE) +
  scale_x_continuous(breaks = 1:3,
                     labels = c("Bayesian: H0", "Bayesian: Inconclusive", "Bayesian: H1"),
                     position = "top", expand = expansion(add = 0.5)) +
  scale_y_continuous(breaks = 1:2, labels = c("Nonsignificant", "Significant"),
                     expand = expansion(add = 0.5)) +
  coord_fixed() +
  labs(title = "Six-Cell Concordance Matrix",
       subtitle = paste0("(Axis not to scale)"),
       x = NULL, y = NULL) +
  theme_project()

ggsave("outputs/figures/concordance.png", concordance_plot, width = 9, height = 6, dpi = 300)


grid_wide <- bf |>
  filter(prior_label %in% c("narrow", "primary", "wide"), is.finite(log10_bf10)) |>
  distinct(claim_id, prior_label, log10_bf10) |>
  tidyr::pivot_wider(names_from = prior_label, values_from = log10_bf10)

spectrum_df <- cl |>
  filter(is.finite(log10_bf10)) |>
  arrange(log10_bf10) |>
  mutate(claim = short_claim(claim_id), row = row_number()) |>
  left_join(grid_wide, by = "claim_id")
for (col in c("narrow", "primary", "wide")) if (!col %in% names(spectrum_df)) spectrum_df[[col]] <- NA_real_
spectrum_df <- spectrum_df |>
  mutate(primary = coalesce(primary, log10_bf10),
         lo = pmin(narrow, primary, wide, na.rm = TRUE),
         hi = pmax(narrow, primary, wide, na.rm = TRUE))

draw_spectrum <- function(df, title, subtitle, xlim) {
  keep <- jeff_log10 >= xlim[1] & jeff_log10 <= xlim[2]
  vl <- jeff_lines[jeff_lines >= xlim[1] & jeff_lines <= xlim[2]]
  ggplot(df) +
    geom_rect(data = jeff_bands,
              aes(xmin = tr_asinh(from), xmax = tr_asinh(to), ymin = -Inf, ymax = Inf, fill = band),
              inherit.aes = FALSE, alpha = 0.50) +
    geom_vline(xintercept = tr_asinh(vl), colour = "white", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
    geom_segment(aes(x = tr_asinh(lo), xend = tr_asinh(hi), y = row, yend = row, colour = concordance_status),
                 linewidth = 1.0, alpha = 0.8) +
    geom_point(aes(tr_asinh(narrow), row), colour = prior_cols[["narrow"]], size = 1.9, na.rm = TRUE) +
    geom_point(aes(tr_asinh(wide), row), colour = prior_cols[["wide"]], size = 1.9, na.rm = TRUE) +
    geom_point(aes(tr_asinh(primary), row), colour = prior_cols[["primary"]], size = 2.8, na.rm = TRUE) +
    scale_colour_manual(values = status_cols, name = "Concordance (segment)", drop = FALSE) +
    scale_fill_manual(values = band_fill, name = "Jeffreys grade") +
    scale_x_continuous(breaks = tr_asinh(jeff_log10[keep]), labels = jeff_labels[keep],
                       sec.axis = sec_axis(~ ., breaks = tr_asinh(jeff_log10[keep]),
                                           labels = sprintf("%.2f", jeff_log10[keep]),
                                           name = expression(log[10](BF[10])))) +
    scale_y_continuous(breaks = df$row, labels = df$claim, expand = expansion(add = 0.6)) +
    coord_cartesian(xlim = tr_asinh(xlim)) +
    guides(fill = guide_legend(nrow = 2, order = 1),
           colour = guide_legend(nrow = 1, order = 2)) +
    labs(title = title, subtitle = subtitle, x = expression(BF[10]~"(Jeffreys scale)"), y = NULL) +
    theme_project() +
    theme(axis.text.y = element_text(size = 7), panel.grid.major.x = element_blank())
}

full_xlim <- range(c(spectrum_df$lo, spectrum_df$hi), na.rm = TRUE) + c(-0.3, 0.3)
evidence_spectrum <- draw_spectrum(
  spectrum_df, "The Evidence Spectrum",
  "Balls = narrow / primary / wide prior scales; segment colour = concordance; shaded band = Jeffreys grade",
  full_xlim)
ggsave("outputs/figures/evidence_spectrum.png", evidence_spectrum,
       width = 10, height = 4 + 0.20 * nrow(spectrum_df), dpi = 300, limitsize = FALSE)

spectrum_zoom_df <- spectrum_df |>
  filter(abs(log10_bf10) <= zoom_lim) |>
  arrange(log10_bf10) |>
  mutate(row = row_number())

evidence_spectrum_zoom <- draw_spectrum(
  spectrum_zoom_df, "The Evidence Spectrum — decision zone",
  paste0("Claims with |log10 BF10| <= ", zoom_lim, "; balls = narrow/primary/wide prior scales, segment colour = concordance"),
  c(-zoom_lim - 0.35, zoom_lim + 0.35))
ggsave("outputs/figures/evidence_spectrum_zoom.png", evidence_spectrum_zoom,
       width = 9, height = 4 + 0.32 * nrow(spectrum_zoom_df), dpi = 300, limitsize = FALSE)

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
  "Axis capped at ±6")
ggsave("outputs/figures/prior_robustness.png", prior_robustness,
       width = 10, height = 4 + 0.22 * nrow(full_df), dpi = 300, limitsize = FALSE)

zoom_df <- prior_wide |>
  dplyr::filter(abs(primary) <= zoom_lim) |>
  dplyr::transmute(claim_id, connector, primary,
                   lo = lower, hi = upper, nar = narrow, pri = primary, wid = wide)
zoom_xlim <- range(c(zoom_df$lo, zoom_df$hi), na.rm = TRUE) + c(-0.1, 0.1)

prior_robustness_zoom <- draw_dumbbell(zoom_df, zoom_xlim,
  "Prior Robustness",
  "Augmented Version")
ggsave("outputs/figures/prior_robustness_zoom.png", prior_robustness_zoom,
       width = 10, height = 4 + 0.30 * nrow(zoom_df), dpi = 300, limitsize = FALSE)


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

slope_labels <- prior_long |>
  group_by(claim_id) |>
  slice_max(order_by = as.integer(prior_label), n = 1, with_ties = FALSE) |>
  ungroup()

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
  ggrepel::geom_text_repel(data = slope_labels, aes(label = short_claim(claim_id)),
                           size = 2.4, direction = "y", hjust = 0, nudge_x = 0.18,
                           segment.size = 0.2, segment.colour = "grey70",
                           max.overlaps = Inf, show.legend = FALSE) +
  coord_cartesian(ylim = c(-6.4, 6.4), clip = "off") +
  labs(title = "Prior Slope Graph",
       subtitle = "Movement Across the Prior Grid; Axis Capped at ±6",
       x = NULL, y = expression(log[10](BF[10]))) +
  theme_project()

ggsave("outputs/figures/prior_slope_graph.png", prior_slope_graph, width = 9, height = 7, dpi = 300)

slope_zoom_lim <- 1.0
prim_lookup <- prior_long |> filter(prior_label == "primary") |> distinct(claim_id, prim = log10_bf10)
slope_zoom_df <- prior_long |> left_join(prim_lookup, by = "claim_id") |> filter(abs(prim) <= slope_zoom_lim)
slope_ylim <- range(slope_zoom_df$log10_bf10, na.rm = TRUE) + c(-0.1, 0.1)
slope_zoom_labels <- slope_zoom_df |>
  group_by(claim_id) |>
  slice_max(order_by = as.integer(prior_label), n = 1, with_ties = FALSE) |>
  ungroup()

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
  ggrepel::geom_text_repel(data = slope_zoom_labels, aes(label = short_claim(claim_id)),
                           size = 2.4, direction = "y", hjust = 0, nudge_x = 0.18,
                           segment.size = 0.2, segment.colour = "grey70",
                           max.overlaps = Inf, show.legend = FALSE) +
  coord_cartesian(ylim = slope_ylim, clip = "off") +
  labs(
    title = "Prior Slope Graph — threshold zoom",
    subtitle = expression("Claims with " * abs(log[10](BF[10])) <= 2),
    x = NULL,
    y = expression(log[10](BF[10]))) +
  theme_project()

ggsave("outputs/figures/prior_slope_graph_zoom.png", prior_slope_graph_zoom, width = 9, height = 6, dpi = 300)


n_lookup <- bf |>
  transmute(claim_id, n_total = safe_numeric(n_total)) |>
  filter(!is.na(n_total)) |>
  distinct(claim_id, .keep_all = TRUE)

fam_df <- cl |>
  filter(is.finite(log10_bf10), !is.na(family)) |>
  left_join(n_lookup, by = "claim_id")

fam_order <- fam_df |>
  group_by(family) |>
  summarise(m = stats::median(log10_bf10), .groups = "drop") |>
  arrange(m) |>
  pull(family)
fam_df <- fam_df |> mutate(family = factor(family, levels = fam_order))

fam_xlim <- range(fam_df$log10_bf10) + c(-0.3, 0.3)
keepf <- jeff_log10 >= fam_xlim[1] & jeff_log10 <= fam_xlim[2]
vlf <- jeff_lines[jeff_lines >= fam_xlim[1] & jeff_lines <= fam_xlim[2]]

fam_layers <- function(with_size) {
  jit <- if (with_size) {
    geom_jitter(aes(size = n_total, colour = concordance_status), height = 0.18, width = 0, alpha = 0.85)
  } else {
    geom_jitter(aes(colour = concordance_status), height = 0.18, width = 0, alpha = 0.85, size = 2.6)
  }
  sz <- if (with_size) scale_size_area(name = "Sample size (N)", max_size = 8) else NULL
  sub <- if (with_size) {
    "Each point a claim on the Jeffreys scale; point size = sample size (N); colour = concordance"
  } else {
    "Each point a claim on the Jeffreys scale; colour = concordance"
  }
  ggplot(fam_df, aes(tr_asinh(log10_bf10), family)) +
    geom_rect(data = jeff_bands, aes(xmin = tr_asinh(from), xmax = tr_asinh(to), ymin = -Inf, ymax = Inf, fill = band),
              inherit.aes = FALSE, alpha = 0.45) +
    geom_vline(xintercept = tr_asinh(vlf), colour = "white", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
    jit +
    scale_fill_manual(values = band_fill, name = "Jeffreys grade") +
    scale_colour_manual(values = status_cols, name = "Concordance", drop = FALSE) +
    sz +
    scale_x_continuous(breaks = tr_asinh(jeff_log10[keepf]), labels = jeff_labels[keepf],
                       sec.axis = sec_axis(~ ., breaks = tr_asinh(jeff_log10[keepf]),
                                           labels = sprintf("%.2f", jeff_log10[keepf]),
                                           name = expression(log[10](BF[10])))) +
    coord_cartesian(xlim = tr_asinh(fam_xlim)) +
    guides(fill = guide_legend(nrow = 2, order = 1),
           colour = guide_legend(nrow = 1, order = 2),
           size = guide_legend(order = 3)) +
    labs(title = "Evidence by Statistical Family", subtitle = sub,
         x = expression(BF[10]~"(Jeffreys scale)"), y = NULL) +
    theme_project() +
    theme(axis.text.y = element_text(size = 8))
}

fam_h <- max(5, 3 + 0.7 * nlevels(fam_df$family))

set.seed(7)
diverging_jeffreys_by_family <- fam_layers(TRUE)
ggsave("outputs/figures/diverging_jeffreys_grades_by_family.png", diverging_jeffreys_by_family,
       width = 10, height = fam_h, dpi = 300, limitsize = FALSE)

set.seed(7)
diverging_jeffreys_by_family_plain <- fam_layers(FALSE)
ggsave("outputs/figures/diverging_jeffreys_grades_by_family_plain.png", diverging_jeffreys_by_family_plain,
       width = 10, height = fam_h, dpi = 300, limitsize = FALSE)


plot_bf_within_frequentist <- function(in_path = "outputs/tables/concordance_claim_level.csv",
                                       out_path = "outputs/figures/bf_within_frequentist.png") {
  d <- readr::read_csv(in_path, show_col_types = FALSE) |>
    dplyr::transmute(panel = dplyr::if_else(frequentist_result == "Significant",
                                            "Significant","Non-significant"),
      verdict = factor(as.character(bf_conclusion),levels = c("H1", "Inconclusive", "H0"))
    ) |>
    dplyr::mutate(panel = factor(panel,levels = c("Significant", "Non-significant"))
    ) |>
    dplyr::count(panel, verdict, .drop = FALSE) |>
    dplyr::group_by(panel) |>
    dplyr::mutate(N = sum(n), proportion = n / N,
                  label = dplyr::if_else(n == 0,"",paste0(n,"\n",scales::percent(proportion, accuracy = 1)))
    ) |> 
    dplyr::ungroup()
  totals <- d |> 
    dplyr::distinct(panel, N)
  cols <- c(H1 = "#08306B", Inconclusive = "#A6A6A6", H0 = "#92C5DE")
  p <- ggplot(d,aes(x = 2,y = proportion,fill = verdict)) +
    geom_col(width = 0.72,colour = "white",linewidth = 1) +
    geom_text(aes(label = label),position = position_stack(vjust = 0.5),colour = "white",fontface = "bold",size = 4,lineheight = 0.95) +
    geom_text(data = totals,aes(x = 0,y = 0,label = paste0("n = ", N)),
              inherit.aes = FALSE,colour = "grey25",fontface = "bold",size = 4.2) +
    facet_wrap(~ panel,nrow = 1) +
    scale_fill_manual(values = cols, breaks = c("H1", "Inconclusive", "H0"),
                      labels = c(H1 = "Evidence for H1: BF10 ≥ 3",Inconclusive = "Inconclusive: 1/3 < BF10 < 3",H0 = "Evidence for H0: BF10 ≤ 1/3"),
                      drop = FALSE,
                      name = NULL) +
    scale_x_continuous(limits = c(0, 2.45),expand = c(0, 0)) +
    coord_polar(theta = "y",start = -pi / 2,clip = "off") +
    labs(title = "Bayesian evidence within frequentist outcomes",subtitle = "") +
    guides(fill = guide_legend(nrow = 1,byrow = TRUE)) +
    theme_void(base_size = 12) +
    theme(plot.background = element_rect(fill = "white",colour = NA),
          panel.background = element_rect(fill = "white",colour = NA),
          legend.background = element_rect(fill = "white",colour = NA),
          legend.key = element_rect(fill = "white",colour = NA),
          plot.title = element_text(face = "bold",size = 17,hjust = 0.5,margin = margin(b = 6)),
          plot.subtitle = element_text(colour = "grey35", size = 10.5,hjust = 0.5,margin = margin(b = 16)),
          strip.text = element_text(face = "bold",size = 13,margin = margin(b = 8)),
          legend.position = "bottom",
          legend.text = element_text(size = 9.5),
          legend.spacing.x = unit(8, "pt"),
          plot.margin = margin(15, 20, 10, 20))
  ggsave(filename = out_path,plot = p,width = 10,height = 5.8,dpi = 300,bg = "white")
  invisible(p)
}
plot_bf_within_frequentist()

plot_bayes_result_proportions <- function(in_path = "outputs/tables/concordance_claim_level.csv",
                                          out_path = "outputs/figures/bayes_result_proportions.png") {
  d <- readr::read_csv(in_path, show_col_types = FALSE) |>
    dplyr::transmute(verdict = factor(as.character(bf_conclusion), levels = c("H1", "Inconclusive", "H0"))) |>
    dplyr::count(verdict, .drop = FALSE) |>
    dplyr::mutate(proportion = n / sum(n),
                  label = dplyr::if_else(n > 0, paste0(n, "\n(", scales::percent(proportion, accuracy = 1), ")"), ""))
  total_n <- sum(d$n)
  cols <- c(H1 = "#08306B", Inconclusive = "#A6A6A6", H0 = "#92C5DE")
  p <- ggplot2::ggplot(d, ggplot2::aes(x = 2, y = proportion, fill = verdict)) +
    ggplot2::geom_col(width = 0.72, colour = "white", linewidth = 1) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      position = ggplot2::position_stack(vjust = 0.5),
      colour = "white", fontface = "bold", size = 4, lineheight = 0.95
    ) +
    ggplot2::annotate(
      "text", x = 0, y = 0, label = paste0("n = ", total_n),
      colour = "grey20", fontface = "bold", size = 4.8
    ) +
    ggplot2::coord_polar(theta = "y", start = -pi / 2) +
    ggplot2::scale_x_continuous(limits = c(0, 2.5), expand = c(0, 0)) +
    scale_fill_manual(values = cols, breaks = c("H1", "Inconclusive", "H0"),
                      labels = c(H1 = "Evidence for H1: BF10 ≥ 3",Inconclusive = "Inconclusive: 1/3 < BF10 < 3",H0 = "Evidence for H0: BF10 ≤ 1/3"),
                      drop = FALSE,
                      name = NULL) +
    ggplot2::labs(
      title = "Distribution of Bayesian conclusions",
      subtitle = ""
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey35", size = 10.5),
      legend.position = "bottom", legend.text = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(15, 20, 15, 20)
    )
  ggplot2::ggsave(out_path, p, width = 7.5, height = 6, dpi = 300, bg = "white")
  invisible(p)
}
plot_bayes_result_proportions()



invisible(list(
  evidence_plane = evidence_plane,
  concordance = concordance_plot,
  evidence_spectrum = evidence_spectrum,
  evidence_spectrum_zoom = evidence_spectrum_zoom,
  prior_robustness = prior_robustness,
  prior_robustness_zoom = prior_robustness_zoom,
  prior_slope_graph = prior_slope_graph,
  prior_slope_graph_zoom = prior_slope_graph_zoom,
  diverging_jeffreys_grades_by_family = diverging_jeffreys_by_family,
  diverging_jeffreys_grades_by_family_plain = diverging_jeffreys_by_family_plain
))


archive_superseded_figures <- function(dir = "outputs/figures", old = "outputs/figures/old",
                                       keep = character()) {
  dir.create(old, showWarnings = FALSE, recursive = TRUE)
  figs <- list.files(dir, pattern = "[.](png|pdf|svg)$")
  move <- setdiff(figs, keep)
  if (length(move) > 0) file.rename(file.path(dir, move), file.path(old, move))
  invisible(move)
}

canonical_figures <- c(
  "evidence_plane.png", "concordance.png",
  "evidence_spectrum.png", "evidence_spectrum_zoom.png",
  "prior_robustness.png", "prior_robustness_zoom.png",
  "prior_slope_graph.png", "prior_slope_graph_zoom.png",
  "diverging_jeffreys_grades_by_family.png",
  "diverging_jeffreys_grades_by_family_plain.png",
  "detailed_evidence_rank_zoom.png", "problematic_claims_pvalues.png",
  "bf_within_frequentist.png", "bayes_result_proportions.png",
  "stability_forest.png",
  "dataset_family_counts.png", "dataset_status_by_family.png",
  "dataset_role_by_family.png", "dataset_sample_size_by_family.png"
)

archive_superseded_figures(keep = canonical_figures)

build_all_figures <- function() {
  invisible(NULL)
}