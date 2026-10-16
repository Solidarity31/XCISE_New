#!/usr/bin/env Rscript
# Figures for the EXCISE 1.5 Python-port benchmark and XCI-skew simulation.
#
# Usage (from the repository root):
#   Rscript simulations/make_figures.R
#
# Inputs  (written by run_speed_benchmark.py / run_xci_skew_sweep.py):
#   simulations/results/speed_benchmark.tsv
#   simulations/results/skew_replicates.tsv
#   simulations/results/skew_snv_level.tsv.gz
# Outputs:
#   simulations/figures/fig{1..5}_*.pdf / .png (300 dpi)
#   simulations/results/fig{1..5}_*_data.csv   (exact numbers plotted)
#   simulations/results/skew_summary.tsv       (mean / SD / range per condition)
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
  optimiser = "#6250d6",  # Python optimiser-only time
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

## ---- Figure 1: runtime Perl vs Python ----------------------------------------
speed <- read_tsv(file.path(res_dir, "speed_benchmark.tsv"), show_col_types = FALSE)
sp_long <- speed %>%
  select(n_cells, seed, n_reads, perl_s, python_s, python_optimiser_s) %>%
  pivot_longer(c(perl_s, python_s, python_optimiser_s), names_to = "what", values_to = "seconds") %>%
  mutate(what = factor(what, c("perl_s", "python_s", "python_optimiser_s"),
                       c("Perl excise1.5.pl (full run)", "Python excise15.py (full run)",
                         "Python optimiser only")))
sp_mean <- sp_long %>% group_by(n_cells, what) %>%
  summarise(seconds = exp(mean(log(pmax(seconds, 1e-3)))), .groups = "drop")  # geometric mean
sp_fold <- speed %>% group_by(n_cells) %>%
  summarise(n_rep = n(), fold = exp(mean(log(speedup))), fold_min = min(speedup),
            fold_max = max(speedup), perl = exp(mean(log(perl_s))),
            identical = sum(bc2xci_identical), .groups = "drop")
write_csv(speed, file.path(res_dir, "fig1_speed_data.csv"))

f1_cols <- c(col[["perl"]], col[["python"]], col[["optimiser"]])
p1a <- ggplot(sp_long, aes(n_cells, seconds, colour = what, shape = what)) +
  geom_line(data = sp_mean, aes(linetype = what), linewidth = 0.5) +
  geom_point(size = 1.6, alpha = 0.85, position = position_jitter(width = 0.02, height = 0, seed = 1)) +
  geom_text(data = sp_fold, inherit.aes = FALSE,
            aes(n_cells, perl * 1.9, label = sprintf("%.1f×", fold)),
            size = 2.5, colour = "black") +
  scale_x_log10(breaks = sort(unique(speed$n_cells)), labels = comma) +
  scale_y_log10(labels = label_number(accuracy = 0.01, drop0trailing = TRUE),
                expand = expansion(mult = c(0.05, 0.15))) +
  scale_colour_manual(values = f1_cols) +
  scale_shape_manual(values = c(16, 17, 1)) +
  scale_linetype_manual(values = c("solid", "solid", "22")) +
  labs(x = "Cells (2,000 SNVs; ~60 reads per cell)", y = "Wall-clock time (s, log scale)",
       title = "Runtime", subtitle = "Numbers: mean speed-up of the full Python run over Perl") +
  theme(legend.position = "bottom", legend.direction = "vertical")

p1b <- ggplot(speed, aes(n_cells, speedup)) +
  geom_hline(yintercept = 1, linewidth = 0.3, colour = col[["ink"]], linetype = "22") +
  geom_line(data = sp_fold, aes(y = fold), colour = col[["python"]], linewidth = 0.5) +
  geom_point(colour = col[["python"]], size = 1.6, alpha = 0.85) +
  geom_text(data = sp_fold, aes(y = fold_max, label = paste0("n=", n_rep)),
            vjust = -1, size = 2.3, colour = col[["ink"]]) +
  scale_x_log10(breaks = sort(unique(speed$n_cells)), labels = comma) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Cells", y = "Speed-up (Perl time / Python time)", title = "Speed-up",
       subtitle = "Points: replicate datasets; line: geometric mean")

fig1 <- (p1a | p1b) + plot_layout(widths = c(1.25, 1)) +
  plot_annotation(
    tag_levels = "A",
    caption = sprintf(paste0(
      "Both tools: -t 4 -j 4 -pretrain, run sequentially on the same simulated BAM. Python time includes BAM reading ",
      "and one-off numba compilation.\nbc2xci output byte-identical between Perl and Python in %d of %d runs."),
      sum(speed$bc2xci_identical), nrow(speed)))
save_fig(fig1, "fig1_runtime_perl_vs_python", 180, 95)

## ---- skew sweep data ---------------------------------------------------------
rep <- read_tsv(file.path(res_dir, "skew_replicates.tsv"), show_col_types = FALSE) %>%
  mutate(ratio = factor(ratio, ratio_levels))
n_rep <- rep %>% count(ratio)
n_lab <- sprintf("%d replicates per ratio; %d cells (5%% Both), %d SNVs",
                 max(n_rep$n), rep$n_cells[1], 1000L)

metrics <- c("hap_acc", "snv_unassigned", "snv_unassigned_nonesc", "escape_unassigned",
             "cell_acc", "cell_acc_informative", "minority_recall", "minority_precision",
             "majority_recall", "recall_Both", "macro_recall", "informative_recovery",
             "mono_recovery", "est_minority_share_of_mono", "wall_s", "optimiser_s")
summary_tbl <- bind_rows(lapply(metrics, function(m) summ(rep, m) %>% mutate(metric = m))) %>%
  select(metric, ratio, n_rep, mean, sd, min, max) %>% arrange(metric, ratio)
write_tsv(summary_tbl, file.path(res_dir, "skew_summary.tsv"))

# replicate points + mean +/- SD, one or two metrics per panel
rep_panel <- function(df, vars, labels, ylab, title, ylim = NULL) {
  d <- df %>% select(ratio, rep, all_of(vars)) %>%
    pivot_longer(all_of(vars), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, vars, labels))
  s <- d %>% group_by(ratio, metric) %>%
    summarise(mean = mean(value), sd = sd(value), .groups = "drop")
  dodge <- if (length(vars) > 1) 0.45 else 0
  pal <- c(col[["python"]], col[["secondary"]])[seq_along(vars)]
  p <- ggplot(d, aes(ratio, value, colour = metric, shape = metric)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = dodge,
                                               jitter.height = 0, seed = 2),
               size = 0.9, alpha = 0.45) +
    geom_errorbar(data = s, inherit.aes = FALSE,
                  aes(ratio, ymin = mean - sd, ymax = mean + sd, colour = metric),
                  width = 0.18, linewidth = 0.4, position = position_dodge(width = dodge)) +
    geom_point(data = s, aes(y = mean), size = 2, position = position_dodge(width = dodge)) +
    scale_colour_manual(values = pal) +
    scale_shape_manual(values = c(16, 17)[seq_along(vars)]) +
    labs(x = "True X1:X2 ratio among monoallelic cells", y = ylab, title = title)
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  if (length(vars) == 1) p <- p + theme(legend.position = "none")
  else p <- p + theme(legend.position = "bottom")
  list(plot = p, data = d)
}

## ---- Figure 2: haplotype accuracy ---------------------------------------------
f2a <- rep_panel(rep, "hap_acc", "Haplotype accuracy",
                 "Haplotype accuracy\n(phased non-escape SNVs, label-switch corrected)",
                 "Haplotype accuracy")
f2b <- rep_panel(rep, c("snv_unassigned_nonesc", "escape_unassigned"),
                 c("Non-escape SNVs left unassigned", "Escape SNVs left unassigned (expected)"),
                 "Fraction of SNVs (X1A = Unk)", "Unassigned SNVs")
write_csv(bind_rows(f2a$data, f2b$data), file.path(res_dir, "fig2_haplotype_data.csv"))
fig2 <- (f2a$plot | f2b$plot) +
  plot_annotation(tag_levels = "A", caption = paste0(
    n_lab, ". Small points: replicates; large points and bars: mean ± SD.\n",
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
    n_lab, ". Small points: replicates; large points and bars: mean ± SD.\n",
    "B: n = true minority (X2) cells per replicate; Low_coverage/Unknown calls count as errors. ",
    "Labels mapped to truth with the global orientation from Fig. 2."))
save_fig(fig3, "fig3_cell_accuracy_vs_skew", 180, 92)

## ---- Figure 4: informative cell recovery ---------------------------------------
f4 <- rep_panel(rep, c("informative_recovery", "mono_recovery"),
                c("All cells called X1, X2 or Both", "True X1/X2 cells called X1 or X2"),
                "Proportion of cells recovered", "Informative cell recovery")
write_csv(f4$data, file.path(res_dir, "fig4_informative_recovery_data.csv"))
fig4 <- f4$plot + labs(caption = paste0(n_lab, ".\nSmall points: replicates; large points and bars: mean ± SD."))
save_fig(fig4, "fig4_informative_cell_recovery", 110, 90)

## ---- Figure 5: haplotypes along chrX ------------------------------------------
snv <- read_tsv(file.path(res_dir, "skew_snv_level.tsv.gz"), show_col_types = FALSE,
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
  geom_point(data = issues, aes(pos / 1e6, 0, shape = issue, fill = issue),
             colour = col[["incorrect"]], size = 1.6, stroke = 0.5) +
  geom_text(data = counts, aes(x = 0, y = 3.05, label = lab), hjust = 0, size = 2.3,
            colour = col[["ink"]]) +
  facet_wrap(~ratio, ncol = 1, labeller = labeller(ratio = function(x) paste0("X1:X2 = ", x, " (replicate 1)"))) +
  scale_colour_manual(values = allele_cols, drop = FALSE) +
  scale_shape_manual(values = issue_shapes, drop = FALSE) +
  scale_fill_manual(values = issue_cols, drop = FALSE, guide = "none") +
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
         shape = guide_legend(nrow = 2, order = 2)) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank())
save_fig(p5, "fig5_haplotype_tracks_chrX", 180, 110)
