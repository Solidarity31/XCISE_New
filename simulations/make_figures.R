#!/usr/bin/env Rscript
# Figures for the EXCISE 1.5 Python-port benchmark and XCI-skew simulation.
#
# Usage (from the repository root):
#   Rscript simulations/make_figures.R
#
# Inputs  (written by run_speed_benchmark.py / run_xci_skew_sweep.py):
#   simulations/results/speed_benchmark.tsv          (Perl / Python / C++ runtimes)
#   simulations/results/skew_indep_replicates.tsv   (independent seeds; Figs 2-5)
#   simulations/results/skew_indep_snv_level.tsv.gz
#   simulations/results/skew_replicates.tsv         (paired control; summary only)
# Outputs:
#   simulations/figures/fig{1..5}_*.pdf / .png (300 dpi)
#   simulations/results/fig{1..5}_*_data.csv   (exact numbers plotted)
#   simulations/results/skew_indep_summary.tsv, skew_paired_summary.tsv
#                                              (mean / SD / range per condition)
#
# Every plotted value is read from the result files; nothing is hard-coded.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(scales)
})
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
sim_dir <- if (length(script)) dirname(normalizePath(script)) else "simulations"
res_dir <- file.path(sim_dir, "results")
fig_dir <- file.path(sim_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- shared style -----------------------------------------------------------
# Validated categorical slots (light surface); the same meaning keeps the same
# colour in every figure. All skew results come from the Python port, so they
# use the Python colour.
col <- c(
  perl      = "#eb6834",  # Perl excise1.5.pl
  python    = "#2a78d6",  # Python excise15.py (and all sweep results)
  secondary = "#1baf7a",  # second metric within a panel (shape-coded too)
  cpp       = "#1baf7a",  # C++ cpp/excise15 (Fig. 1 only; direct-labelled)
  ref       = "#6250d6",  # REF allele on X1
  alt       = "#eda100",  # ALT allele on X1
  incorrect = "#e34948",
  unknown   = "#8a8984",
  ink       = "#52514e"
)

theme_excise <- function(base = 8) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = col[["ink"]]),
      axis.ticks = element_line(linewidth = 0.3, colour = col[["ink"]]),
      axis.text = element_text(colour = "black"),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey90"),
      legend.key.size = unit(3.5, "mm"),
      legend.title = element_blank(),
      legend.background = element_blank(),
      plot.title = element_text(face = "bold", size = base + 1),
      plot.subtitle = element_text(colour = col[["ink"]], size = base - 0.5),
      plot.caption = element_text(colour = col[["ink"]], size = base - 1.5, hjust = 0),
      plot.tag = element_text(face = "bold", size = base + 2),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base)
    )
}
theme_set(theme_excise())

save_fig <- function(p, name, w_mm, h_mm) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w_mm, height = h_mm,
         units = "mm", device = pdf, useDingbats = FALSE)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w_mm, height = h_mm,
         units = "mm", dpi = 300, bg = "white")
  message("wrote ", file.path(fig_dir, name), ".{pdf,png}")
}

ratio_levels <- c("50:50", "55:45", "70:30", "90:10", "95:5", "99:1")

# mean +/- SD summary used for every replicate panel
summ <- function(df, value) {
  df %>% group_by(ratio) %>%
    summarise(n_rep = n(), mean = mean(.data[[value]]), sd = sd(.data[[value]]),
              min = min(.data[[value]]), max = max(.data[[value]]), .groups = "drop")
}

## ---- Figure 1: runtime, original XCISE vs EXCISE 1.5 (Perl / Python / C++) ----
speed <- read_tsv(file.path(res_dir, "speed_benchmark.tsv"), show_col_types = FALSE)
has_cpp <- "cpp_s" %in% names(speed) && !anyNA(speed$cpp_s)
has_orig <- "xcise_s" %in% names(speed) && !anyNA(speed$xcise_s)
tools <- c(xcise_s = "Original XCISE xcise.pl", perl_s = "EXCISE 1.5 Perl", python_s = "EXCISE 1.5 Python",
           cpp_s = "EXCISE 1.5 C++")
tool_cols <- c(col[["ink"]], col[["perl"]], col[["python"]], col[["cpp"]])
tool_shapes <- c(18, 16, 17, 15)
tool_lty <- c("22", "solid", "solid", "solid")
keep <- c(has_orig, TRUE, TRUE, has_cpp)
tools <- tools[keep]
tool_cols <- setNames(tool_cols[keep], tools); tool_shapes <- setNames(tool_shapes[keep], tools)
tool_lty <- setNames(tool_lty[keep], tools)
base_col <- if (has_orig) "xcise_s" else "perl_s"           # speed-up baseline
base_lab <- tools[[base_col]]
sp_long <- speed %>%
  select(n_cells, seed, all_of(names(tools))) %>%
  pivot_longer(all_of(names(tools)), names_to = "key", values_to = "seconds") %>%
  mutate(what = factor(tools[key], tools))
sp_mean <- sp_long %>% group_by(n_cells, what) %>%
  summarise(seconds = exp(mean(log(seconds))), .groups = "drop")      # geometric mean
fold_long <- sp_long %>%
  left_join(speed %>% select(n_cells, seed, base = all_of(base_col)), by = c("n_cells", "seed")) %>%
  filter(key != base_col) %>%
  mutate(fold = base / seconds, what = droplevels(what))
fold_mean <- fold_long %>% group_by(n_cells, what) %>%
  summarise(n_rep = n(), fold = exp(mean(log(fold))), .groups = "drop")
write_csv(speed, file.path(res_dir, "fig1_speed_data.csv"))
write_csv(fold_mean, file.path(res_dir, "fig1_speedup_geomean.csv"))

xmax <- max(speed$n_cells)
end_lab <- function(d) d %>% filter(n_cells == xmax) %>% mutate(lab = as.character(what))
p1a <- ggplot(sp_long, aes(n_cells, seconds, colour = what, shape = what)) +
  geom_line(data = sp_mean, aes(linetype = what), linewidth = 0.5) +
  geom_point(size = 1.6, alpha = 0.85, position = position_jitter(width = 0.02, height = 0, seed = 1)) +
  geom_text(data = end_lab(sp_mean), aes(label = lab), hjust = 0, nudge_x = 0.06, size = 2.3,
            show.legend = FALSE) +
  scale_x_log10(breaks = sort(unique(speed$n_cells)), labels = comma,
                expand = expansion(mult = c(0.04, 0.55))) +
  scale_y_log10(labels = label_number(accuracy = 0.1, drop0trailing = TRUE, big.mark = ",")) +
  scale_colour_manual(values = tool_cols) + scale_shape_manual(values = tool_shapes) +
  scale_linetype_manual(values = tool_lty) +
  labs(x = "Cells (2,000 SNVs; ~60 reads per cell)", y = "Wall-clock time per run (s, log scale)",
       title = "Runtime", subtitle = "Points: datasets; lines: geometric mean") +
  theme(legend.position = "none")

lab_above <- tools[["python_s"]]                             # Python labels above, others below
p1b <- ggplot(fold_long, aes(n_cells, fold, colour = what, shape = what)) +
  geom_hline(yintercept = 1, linewidth = 0.3, colour = col[["ink"]], linetype = "22") +
  geom_line(data = fold_mean, linewidth = 0.5) +
  geom_point(size = 1.6, alpha = 0.85) +
  geom_text(data = fold_mean, aes(label = ifelse(fold >= 100, comma(round(fold)), sprintf("%.0f", fold)),
                                  vjust = ifelse(what == lab_above, -1.1, 2.0)),
            size = 2.1, show.legend = FALSE) +
  geom_text(data = end_lab(fold_mean), aes(label = sub("EXCISE 1.5 ", "", lab)), hjust = 0, nudge_x = 0.06,
            size = 2.3, show.legend = FALSE) +
  scale_x_log10(breaks = sort(unique(speed$n_cells)), labels = comma,
                expand = expansion(mult = c(0.06, 0.18))) +
  scale_y_log10(breaks = 10^(0:5), labels = label_comma(), expand = expansion(mult = c(0.08, 0.12))) +
  scale_colour_manual(values = tool_cols) + scale_shape_manual(values = tool_shapes) +
  labs(x = "Cells", y = sprintf("Speed-up over %s (log scale)", sub(" xcise.pl", "", base_lab)),
       title = sprintf("Speed-up over %s", sub(" xcise.pl", "", base_lab)),
       subtitle = sprintf("Numbers: geometric mean (x) of %d datasets per size", max(fold_mean$n_rep))) +
  theme(legend.position = "none")

cap1 <- sprintf(paste0(
  "All tools: 4 tries (-t 4) on the same simulated BAMs (%d datasets), run one at a time; times are full runs ",
  "including BAM reading.\nEXCISE 1.5: -j 4 -pretrain.%s\n",
  "Python includes one-off numba compilation. Optimiser alone: Python %s s, C++ %s s.\n",
  "Cell accuracy vs simulated truth (mean): %s.\n",
  "bc2xci identical to Perl EXCISE 1.5: Python %d/%d, C++ %d/%d runs; the rest reached different local optima ",
  "(at most %d cell labelled differently)."),
  nrow(speed),
  if (has_orig) " Original XCISE (Vityay/XCISE @3a1d063): -u 4 (same SNV filter); no -j or -pretrain, tries run one after another." else "",
  paste(sprintf("%.2f", range(speed$python_optimiser_s)), collapse = "-"),
  if (has_cpp) sprintf("at most %.2f", max(speed$cpp_optimiser_s)) else "NA",
  paste(sprintf("%s %.3f", sub("EXCISE 1.5 |Original XCISE xcise.pl", "", c(if (has_orig) "Original" , "Perl", "Python", if (has_cpp) "C++")),
                c(if (has_orig) mean(speed$xcise_cell_acc), mean(speed$perl_cell_acc), mean(speed$python_cell_acc),
                  if (has_cpp) mean(speed$cpp_cell_acc))), collapse = ", "),
  sum(speed$bc2xci_identical), nrow(speed),
  if (has_cpp) sum(speed$cpp_bc2xci_identical_perl) else NA, nrow(speed),
  max(c(speed$cells_label_differs, if (has_cpp) speed$cpp_cells_label_differs_perl)))
fig1 <- (p1a | p1b) + plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = "A", caption = cap1)
save_fig(fig1, "fig1_runtime_benchmark", 180, 108)

## ---- skew sweep data ---------------------------------------------------------
read_rep <- function(tag) read_tsv(file.path(res_dir, paste0(tag, "_replicates.tsv")),
                                   show_col_types = FALSE) %>%
  mutate(ratio = factor(ratio, ratio_levels))
rep <- read_rep("skew_indep")          # main analysis: independent seeds
rep_paired <- read_rep("skew")         # control: common random numbers
n_rep <- rep %>% count(ratio)
n_lab <- sprintf("%d independent replicates per ratio; %d cells (5%% Both), %d SNVs",
                 max(n_rep$n), rep$n_cells[1], 1000L)

metrics <- c("hap_acc", "snv_unassigned", "snv_unassigned_nonesc", "escape_unassigned",
             "cell_acc", "cell_acc_informative", "minority_recall", "minority_precision",
             "majority_recall", "recall_Both", "macro_recall", "informative_recovery",
             "mono_recovery", "est_minority_share_of_mono", "wall_s", "optimiser_s")
summ_all <- function(df) bind_rows(lapply(metrics, function(m) summ(df, m) %>% mutate(metric = m))) %>%
  select(metric, ratio, n_rep, mean, sd, min, max) %>% arrange(metric, ratio)
write_tsv(summ_all(rep), file.path(res_dir, "skew_indep_summary.tsv"))
write_tsv(summ_all(rep_paired), file.path(res_dir, "skew_paired_summary.tsv"))

# replicate points + mean +/- SD, one or two metrics per panel
rep_panel <- function(df, vars, labels, ylab, title, ylim = NULL) {
  d <- df %>% select(ratio, rep, all_of(vars)) %>%
    pivot_longer(all_of(vars), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, vars, labels))
  s <- d %>% group_by(ratio, metric) %>%
    summarise(mean = mean(value), lo = min(value), hi = max(value), .groups = "drop")
  dodge <- if (length(vars) > 1) 0.45 else 0
  pal <- c(col[["python"]], col[["secondary"]])[seq_along(vars)]
  p <- ggplot(d, aes(ratio, value, colour = metric, shape = metric)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = dodge,
                                               jitter.height = 0, seed = 2),
               size = 0.9, alpha = 0.45) +
    geom_errorbar(data = s, inherit.aes = FALSE,
                  aes(ratio, ymin = lo, ymax = hi, colour = metric),
                  width = 0.18, linewidth = 0.4, position = position_dodge(width = dodge)) +
    geom_point(data = s, aes(y = mean), size = 2, position = position_dodge(width = dodge)) +
    scale_colour_manual(values = pal) +
    scale_shape_manual(values = c(16, 17)[seq_along(vars)]) +
    labs(x = "True X1:X2 ratio among monoallelic cells", y = ylab, title = title)
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  if (length(vars) == 1) p <- p + theme(legend.position = "none")
  else p <- p + theme(legend.position = "bottom", legend.direction = "vertical")
  list(plot = p, data = d)
}

## ---- Figure 2: haplotype accuracy ---------------------------------------------
f2a <- rep_panel(rep, "hap_acc", "Haplotype accuracy",
                 "Haplotype accuracy\n(phased non-escape SNVs, label-switch corrected)",
                 "Haplotype accuracy")
f2a$plot <- f2a$plot +
  scale_y_continuous(limits = c(min(0.95, min(rep$hap_acc)), 1), breaks = seq(0.95, 1, 0.01)) +
  labs(subtitle = sprintf("All %d runs: accuracy %.4f-%.4f\n(%s-%s non-escape SNVs phased per run)",
                          nrow(rep), min(rep$hap_acc), max(rep$hap_acc),
                          comma(min(rep$snvs_phased_nonesc)), comma(max(rep$snvs_phased_nonesc))))
rep_pct <- rep %>% mutate(unassigned_nonesc_pct = 100 * snv_unassigned_nonesc)
f2b <- rep_panel(rep_pct, "unassigned_nonesc_pct", "Non-escape SNVs left unassigned",
                 "Non-escape SNVs left unassigned (%)", "Unassigned SNVs")
f2b$plot <- f2b$plot + scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(subtitle = sprintf("Escape SNVs left unassigned (desired outcome):\n%.0f-%.0f%% across runs",
                          100 * min(rep$escape_unassigned), 100 * max(rep$escape_unassigned)))
write_csv(bind_rows(f2a$data, f2b$data), file.path(res_dir, "fig2_haplotype_data.csv"))
fig2 <- (f2a$plot | f2b$plot) +
  plot_annotation(tag_levels = "A", caption = paste0(
    n_lab, ". Small points: replicates; large points and bars: mean and min-max range.\n",
    "Accuracy uses one global X1/X2 orientation per run, chosen from the SNV haplotypes ",
    "(EXCISE labels are arbitrary)."))
save_fig(fig2, "fig2_haplotype_accuracy_vs_skew", 180, 90)

## ---- Figure 3: cell assignment & minority recall -------------------------------
f3a <- rep_panel(rep, c("cell_acc", "macro_recall"),
                 c("Accuracy (all cells)", "Mean per-class recall (X1, X2, Both)"),
                 "Fraction", "Cell assignment")
f3b <- rep_panel(rep, c("minority_recall", "minority_precision"),
                 c("Minority-class recall", "Minority-class precision"),
                 "Fraction", "Minority class (true X2 cells)")
nmin <- rep %>% group_by(ratio) %>%
  summarise(lab = if (n_distinct(n_minority) == 1) sprintf("n=%d", n_minority[1])
                  else sprintf("n=%d-%d", min(n_minority), max(n_minority)), .groups = "drop")
f3b$plot <- f3b$plot +
  geom_text(data = nmin, inherit.aes = FALSE, aes(ratio, Inf, label = lab),
            vjust = 1.2, size = 2.3, colour = col[["ink"]]) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))
write_csv(bind_rows(f3a$data, f3b$data), file.path(res_dir, "fig3_cell_assignment_data.csv"))
fig3 <- (f3a$plot | f3b$plot) +
  plot_annotation(tag_levels = "A", caption = paste0(
    n_lab, ". Small points: replicates; large points and bars: mean and min-max range.\n",
    "B: n = true minority (X2) cells per replicate; Low_coverage/Unknown calls count as errors. ",
    "Labels mapped to truth with the global orientation from Fig. 2."))
save_fig(fig3, "fig3_cell_accuracy_vs_skew", 180, 92)

## ---- Figure 4: informative cell recovery ---------------------------------------
f4 <- rep_panel(rep, c("informative_recovery", "mono_recovery"),
                c("All cells called X1, X2 or Both", "True X1/X2 cells called X1 or X2"),
                "Proportion of cells recovered", "Informative cell recovery")
write_csv(f4$data, file.path(res_dir, "fig4_informative_recovery_data.csv"))
fig4 <- f4$plot + labs(caption = paste0(n_lab, ".\nSmall points: replicates; large points and bars: mean and min-max range."))
save_fig(fig4, "fig4_informative_cell_recovery", 110, 90)

## ---- Figure 5: haplotypes along chrX ------------------------------------------
snv <- read_tsv(file.path(res_dir, "skew_indep_snv_level.tsv.gz"), show_col_types = FALSE,
                col_types = cols(inferred_x1_allele = col_character())) %>%
  mutate(ratio = factor(ratio, ratio_levels))
sel <- snv %>% filter(ratio %in% c("50:50", "99:1"), rep == 1)
tracks <- bind_rows(
  sel %>% transmute(ratio, pos, track = "Truth",
                    allele = ifelse(escape == 1, "Escape (biallelic)",
                                    ifelse(true_x1_allele == 1, "REF on X1", "ALT on X1"))),
  sel %>% transmute(ratio, pos, track = "EXCISE",
                    allele = case_when(status == "filtered" ~ "Filtered (not tested)",
                                       inferred_x1_allele == "1" ~ "REF on X1",
                                       inferred_x1_allele == "2" ~ "ALT on X1",
                                       TRUE ~ "Unassigned (Unk)"))) %>%
  mutate(track = factor(track, c("EXCISE", "Truth")),
         allele = factor(allele, c("REF on X1", "ALT on X1", "Unassigned (Unk)",
                                   "Escape (biallelic)", "Filtered (not tested)")))
issues <- sel %>% filter(status %in% c("incorrect", "unassigned", "escape_phased", "filtered")) %>%
  mutate(issue = factor(status, c("incorrect", "unassigned", "escape_phased", "filtered"),
                        c("Wrong orientation", "Unassigned (non-escape)",
                          "Escape SNV phased", "Filtered")))
counts <- sel %>% group_by(ratio) %>%
  summarise(lab = sprintf("correct %d | wrong %d | unassigned %d | escape kept Unk %d/%d",
                          sum(status == "correct"), sum(status == "incorrect"),
                          sum(status == "unassigned"), sum(status == "escape_unassigned"),
                          sum(escape == 1)), .groups = "drop")
write_csv(sel, file.path(res_dir, "fig5_snv_tracks_data.csv"))

allele_cols <- c("REF on X1" = col[["ref"]], "ALT on X1" = col[["alt"]],
                 "Unassigned (Unk)" = col[["unknown"]], "Escape (biallelic)" = "grey80",
                 "Filtered (not tested)" = "grey92")
issue_shapes <- c("Wrong orientation" = 4, "Unassigned (non-escape)" = 1,
                  "Escape SNV phased" = 2, "Filtered" = 3)
issue_cols <- c("Wrong orientation" = col[["incorrect"]], "Unassigned (non-escape)" = col[["unknown"]],
                "Escape SNV phased" = col[["incorrect"]], "Filtered" = col[["unknown"]])

p5 <- ggplot() +
  geom_segment(data = tracks, aes(x = pos / 1e6, xend = pos / 1e6,
                                  y = as.numeric(track) - 0.38, yend = as.numeric(track) + 0.38,
                                  colour = allele), linewidth = 0.18) +
  geom_point(data = issues, aes(pos / 1e6, 0, shape = issue),
             colour = issue_cols[as.character(issues$issue)], size = 1.8, stroke = 0.6) +
  geom_text(data = counts, aes(x = 0, y = 3.05, label = lab), hjust = 0, size = 2.3,
            colour = col[["ink"]]) +
  facet_wrap(~ratio, ncol = 1, labeller = labeller(ratio = function(x) paste0("X1:X2 = ", x, " (replicate 1)"))) +
  scale_colour_manual(values = allele_cols, drop = FALSE) +
  scale_shape_manual(values = issue_shapes, drop = FALSE) +
  scale_y_continuous(breaks = c(0, 1, 2), labels = c("Errors /\nunassigned", "EXCISE", "Truth"),
                     limits = c(-0.6, 3.3), expand = c(0, 0)) +
  scale_x_continuous(expand = expansion(mult = 0.01)) +
  labs(x = "Position on simulated chrX (Mb)", y = NULL,
       title = "Inferred vs true haplotypes along chromosome X",
       caption = paste0("Each tick is one SNV, coloured by which allele lies on haplotype X1 ",
                        "(EXCISE orientation corrected by the global label switch).\n",
                        "Bottom row marks every SNV that is wrongly oriented, left unassigned, ",
                        "or removed by the read-depth filter.")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.5), nrow = 2, order = 1),
         shape = guide_legend(nrow = 2, order = 2,
                              override.aes = list(colour = unname(issue_cols)))) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank())
save_fig(p5, "fig5_haplotype_tracks_chrX", 180, 110)

## ---- Figure 6: parameter sweep (-u x -p) ------------------------------------
ps_path <- file.path(res_dir, "param_sweep_replicates.tsv")
if (file.exists(ps_path)) {
  ps <- read_tsv(ps_path, show_col_types = FALSE)
  regime_lab <- c(standard = "Standard\n40 reads/cell, 2% noise",
                  low_coverage = "Low coverage\n10 reads/cell",
                  high_noise = "High noise\n10% inactive-X reads",
                  sparse_snvs = "Sparse SNVs\n3,000 SNVs, 15 reads/cell")
  n_mono <- ps$n_X1[1] + ps$n_X2[1]
  pss <- ps %>%
    mutate(wrong_x_pct = 100 * wrong_x_cells / (n_X1 + n_X2),
           unassigned_pct = 100 * snv_unassigned_nonesc) %>%
    group_by(regime, min_umis, penalty) %>%
    summarise(n_rep = n(),
              across(c(cell_acc, wrong_x_pct, unassigned_pct, mono_recovery, recall_Both, hap_acc,
                       snvs_pass_filter, escape_unassigned),
                     list(mean = mean, sd = sd, min = min, max = max)), .groups = "drop") %>%
    mutate(regime = factor(regime, names(regime_lab)))
  # tables keep the plain regime keys; the plots use the multi-line labels
  write_tsv(pss %>% arrange(regime, min_umis, penalty), file.path(res_dir, "param_sweep_summary.tsv"))
  write_csv(pss %>% select(regime, min_umis, penalty, n_rep, cell_acc_mean, wrong_x_pct_mean,
                           unassigned_pct_mean), file.path(res_dir, "fig6_param_sweep_data.csv"))
  pss <- pss %>%
    mutate(regime = factor(regime_lab[as.character(regime)], regime_lab),
           u = factor(min_umis, sort(unique(min_umis))), p = factor(penalty, sort(unique(penalty))))
  dflt <- pss %>% filter(min_umis == 4, penalty == 5)
  heat <- function(var, title, fmt, low, high, better_high) {
    ggplot(pss, aes(p, u, fill = .data[[var]])) +
      geom_tile(colour = "white", linewidth = 0.4) +
      geom_tile(data = dflt, fill = NA, colour = "black", linewidth = 0.6) +
      # white text only on the darker half of the colour ramp
      geom_text(aes(label = sprintf(fmt, .data[[var]]),
                    colour = (.data[[var]] - min(.data[[var]])) / diff(range(.data[[var]])) > 0.55),
                size = 1.7, show.legend = FALSE) +
      scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
      scale_fill_gradient(low = low, high = high, name = NULL) +
      facet_wrap(~regime, nrow = 1) +
      labs(x = "Discordant penalty (-p)", y = "Min UMIs per SNV (-u)", title = title) +
      theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
            legend.key.height = unit(6, "mm"), legend.key.width = unit(2.5, "mm"))
  }
  f6a <- heat("cell_acc_mean", "Cell accuracy (all cells; Low_coverage/Unknown count as errors)", "%.2f",
              "#f2f6fc", "#1f5fae", TRUE)
  f6b <- heat("wrong_x_pct_mean", "Cells assigned to the wrong X (% of true X1/X2 cells)", "%.1f",
              "#fdf3f2", "#b42c2b", FALSE)
  f6c <- heat("unassigned_pct_mean", "Non-escape SNVs left unassigned (%)", "%.0f",
              "#f4f4f2", "#52514e", FALSE)
  fig6 <- (f6a / f6b / f6c) +
    plot_annotation(tag_levels = "A", caption = sprintf(paste0(
      "Mean of %d simulated replicates per regime and setting (2,000 cells: 5%% Both, X1:X2 = 50:50; 5%% escape SNVs). ",
      "EXCISE 1.5 (cpp/excise15, identical to excise1.5.pl)\nwith -t 100 -pretrain; ",
      "only -u and -p vary. Black outline: defaults (-u 4, -p 5). Haplotype accuracy and the other metrics are in ",
      "results/param_sweep_summary.tsv."), max(pss$n_rep)))
  save_fig(fig6, "fig6_param_sweep_min_umis_penalty", 180, 215)
}
