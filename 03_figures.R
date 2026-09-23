# Reproduce the four simulation figures in the paper.
# Run from this folder: Rscript 03_figures.R
# Default: read the supplied data/paper summaries; write four PDFs to figures/.
# Optional: --png also writes 320-dpi PNGs; --out-dir=PATH changes the destination.
# Fresh fits: Rscript 03_figures.R --memberships=results/memberships.csv --png
# Recompute archived results: --memberships=data/paper/memberships.csv.gz
# Or pass --input-dir=PATH containing the same three summary CSVs as data/paper/.
# Install ggplot2, dplyr and scales; fresh memberships additionally need mclust.
# Fresh data require all three effects, all three intervals, 60 subjects per fit,
# and >=2 matched replications per effect. Columns: signal,effect,delta,replication,
# subject,truth,ct_group,dt_group. Signals are Weak/Mid/Large; effect=.3/.6/.9;
# delta=.5/1/5; subject=1:60; truth=rep(1:3,each=20). Cluster labels may be arbitrary.
# For a different design, adjust these checks, effect labels and axis scales.
# Keep the original settings below for the paper: strict ARI > .70, 95% Wilson
# intervals for proportions, t-based 95% Monte Carlo CIs for replication means,
# paired full-precision ARIs, and pair classifications matched over all 3 rates.
# Frozen summaries retain the original report's precision. Fonts/rendering can
# vary with R/ggplot2/system versions; PDFs require an R build with Cairo support.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(scales) })
args <- commandArgs(trailingOnly = TRUE)
if (any(!grepl("^(--png|--out-dir=.+|--input-dir=.+|--memberships=.+)$", args)))
  stop("Options: --png --out-dir=PATH --input-dir=PATH --memberships=FILE")
option <- function(name, default) {
  value <- args[startsWith(args, paste0("--", name, "="))]
  if (length(value) > 1L) stop("Repeated option: ", name)
  if (length(value)) sub(paste0("^--", name, "="), "", value) else default
}
input_dir <- option("input-dir", "data/paper")
out_dir <- option("out-dir", "figures")
memberships_file <- option("memberships", "")
if (nzchar(memberships_file) && any(startsWith(args, "--input-dir=")))
  stop("Choose --memberships or --input-dir, not both.")
if (!capabilities("cairo")) stop("PDF export requires Cairo support in R.")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
effect_lookup <- c(Weak = .30, Mid = .60, Large = .90)
effect_labels <- c(Weak = "Weak (0.30)", Mid = "Moderate (0.60)", Large = "Large (0.90)")
method_levels <- c("ctS-GIMME", "DT-GIMME (S-GIMME)")
summary_files <- c("ari_summary.csv", "paired_ari_summary.csv", "cross_rate_stability_summary.csv")

# Recompute just the statistics needed by the active figures from new fits.
mean_ci <- function(x) {
  n <- length(x); se <- sd(x) / sqrt(n); half <- qt(.975, n - 1L) * se
  data.frame(n = n, mean = mean(x), sd = sd(x), mcse = se,
             ci_low = mean(x) - half, ci_high = mean(x) + half)
}
wilson <- function(x) {
  n <- length(x); p <- mean(x); z <- qnorm(.975); denominator <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denominator
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denominator
  data.frame(prop_gt_070 = p, prop_070_ci_low = max(0, center - half),
             prop_070_ci_high = min(1, center + half))
}
summarize_groups <- function(data, keys, fun) {
  bind_rows(lapply(split(data, interaction(data[keys], drop = TRUE)), function(d) {
    values <- fun(d); labels <- d[rep(1L, nrow(values)), keys, drop = FALSE]
    rownames(labels) <- NULL
    cbind(labels, values)
  }))
}
if (nzchar(memberships_file)) {
  if (!requireNamespace("mclust", quietly = TRUE)) stop("Install mclust for fresh memberships.")
  members <- read.csv(memberships_file, check.names = FALSE)
  required <- c("signal", "effect", "delta", "replication", "subject", "truth", "ct_group", "dt_group")
  if (!all(required %in% names(members)) || !nrow(members)) stop("Invalid membership CSV schema.")
  members <- members[required]
  if (anyNA(members) || any(trimws(as.character(members$ct_group)) == "") ||
      any(trimws(as.character(members$dt_group)) == "")) stop("Memberships must be complete.")
  numeric_fields <- c("effect", "delta", "replication", "subject", "truth")
  if (any(!vapply(members[numeric_fields], function(x) is.numeric(x) && all(is.finite(x)), logical(1))))
    stop("Effect, delta, replication, subject and truth must be finite numbers.")
  if (!setequal(unique(members$signal), names(effect_lookup)) ||
      any(members$effect != unname(effect_lookup[members$signal])) ||
      !setequal(unique(members$delta), c(.5, 1, 5)) ||
      any(members$replication < 1 | members$replication != as.integer(members$replication)))
    stop("Expected all three paper effects and intervals, with positive integer replication IDs.")
  if (anyDuplicated(members[c("signal", "delta", "replication", "subject")]))
    stop("Duplicate subject within a condition and replication.")
  fits <- split(members, interaction(members[c("signal", "delta", "replication")], drop = TRUE))
  if (any(!vapply(fits, function(d) nrow(d) == 60L && setequal(d$subject, 1:60) &&
                 all(d$truth == ceiling(d$subject / 20)), logical(1))))
    stop("Every fit must contain subjects 1:60 with three consecutive true groups of 20.")
  blocks <- split(members, interaction(members[c("signal", "replication")], drop = TRUE))
  if (any(!vapply(blocks, function(d) nrow(d) == 180L && setequal(d$delta, c(.5, 1, 5)), logical(1))))
    stop("Every replication must contain all three matched intervals; incomplete fits cannot be plotted.")
  if (any(vapply(split(members$replication, members$signal), function(x) length(unique(x)), integer(1)) < 2L))
    stop("At least two complete replications per effect are needed for confidence intervals.")
  ari <- bind_rows(lapply(fits, function(d) data.frame(
    signal = d$signal[1L], effect = d$effect[1L], delta = d$delta[1L], replication = d$replication[1L],
    ctARI = mclust::adjustedRandIndex(d$ct_group, d$truth),
    dtARI = mclust::adjustedRandIndex(d$dt_group, d$truth))))
  if (any(!is.finite(as.matrix(ari[c("ctARI", "dtARI")])))) stop("Non-finite ARI.")
  ari_long <- bind_rows(lapply(seq_along(method_levels), function(i)
    data.frame(ari[c("signal", "effect", "delta", "replication")],
               method = method_levels[i], ARI = ari[[c("ctARI", "dtARI")[i]]])))
  ari_summary <- summarize_groups(ari_long, c("signal", "effect", "delta", "method"),
                                  function(d) cbind(mean_ci(d$ARI), wilson(d$ARI > .70)))
  paired_ari_summary <- summarize_groups(ari, c("signal", "effect", "delta"),
                                         function(d) mean_ci(d$ctARI - d$dtARI))
  pair_index <- t(combn(1:60, 2L)); truth <- rep(1:3, each = 20L)
  same_truth <- truth[pair_index[, 1L]] == truth[pair_index[, 2L]]
  pair_metrics <- bind_rows(lapply(blocks, function(d) {
    rates <- lapply(c(.5, 1, 5), function(delta) {
      z <- d[d$delta == delta, ]; z[order(z$subject), ]
    })
    bind_rows(lapply(seq_along(method_levels), function(i) {
      count <- rowSums(vapply(rates, function(z) {
        g <- z[[c("ct_group", "dt_group")[i]]]
        g[pair_index[, 1L]] == g[pair_index[, 2L]]
      }, logical(nrow(pair_index))))
      bind_rows(lapply(c(TRUE, FALSE), function(same) data.frame(
        signal = d$signal[1L], effect = d$effect[1L], replication = d$replication[1L],
        method = method_levels[i], relation = if (same) "True same-group pairs" else "True different-group pairs",
        stable_correct = mean(count[same_truth == same] == if (same) 3L else 0L),
        unstable = mean(count[same_truth == same] > 0L & count[same_truth == same] < 3L))))
    }))
  }))
  cross_rate_stability_summary <- summarize_groups(pair_metrics, c("signal", "effect", "method", "relation"),
    function(d) bind_rows(lapply(c("stable_correct", "unstable"), function(metric)
      cbind(metric = metric, mean_ci(d[[metric]])))))
  summaries <- list(ari_summary, paired_ari_summary, cross_rate_stability_summary)
  summaries <- lapply(summaries, function(d) { d$effect_label <- unname(effect_labels[d$signal]); d })
  dir.create(file.path(out_dir, "summaries"), showWarnings = FALSE)
  for (i in seq_along(summaries)) write.csv(summaries[[i]], file.path(out_dir, "summaries", summary_files[i]), row.names = FALSE)
} else summaries <- lapply(file.path(input_dir, summary_files), read.csv, check.names = FALSE)

# Restore the report's factor order after reading CSVs.
summaries <- lapply(summaries, function(d) {
  d$signal <- factor(d$signal, levels = names(effect_lookup))
  d$effect_label <- factor(d$effect_label, levels = unname(effect_labels))
  if ("method" %in% names(d)) d$method <- factor(d$method, levels = method_levels)
  if ("relation" %in% names(d)) d$relation <- factor(d$relation, levels = c("True same-group pairs", "True different-group pairs"))
  if ("metric" %in% names(d)) d$metric <- factor(d$metric, levels = c("stable_correct", "unstable", "stable_wrong"))
  if (anyNA(d$effect_label) || ("method" %in% names(d) && anyNA(d$method))) stop("Unknown effect or method label.")
  d |> arrange(across(any_of(c("signal", "delta", "method", "relation", "metric"))))
})
ari_summary <- summaries[[1L]]; paired_ari_summary <- summaries[[2L]]
cross_rate_stability_summary <- summaries[[3L]]
aggie_blue <- "#022851"; aggie_gold <- "#FFBF00"; neutral_charcoal <- "#595959"
palette_methods <- setNames(c(aggie_blue, aggie_gold), method_levels)
theme_report <- function(base_size = 10.5) {
  theme_classic(base_size = base_size, base_family = "sans") + theme(
    strip.text = element_text(face = "bold", color = "black", margin = margin(4, 4, 4, 4)),
    strip.background = element_blank(), axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"), axis.line = element_line(color = "black", linewidth = .35),
    axis.ticks = element_line(color = "black", linewidth = .35), legend.position = "top",
    legend.justification = "left", legend.box.just = "left", legend.title = element_blank(),
    legend.key = element_blank(), legend.key.width = grid::unit(1.2, "lines"),
    panel.spacing = grid::unit(.9, "lines"), plot.margin = margin(8, 10, 8, 8))
}
method_scales <- function(labels = waiver()) list(
  scale_color_manual(values = palette_methods, labels = labels),
  scale_linetype_manual(values = setNames(c("solid", "22"), method_levels), labels = labels),
  scale_shape_manual(values = setNames(c(16, 17), method_levels), labels = labels))
delta_scale <- function() scale_x_continuous(breaks = c(.5, 1, 5), labels = c("0.5", "1", "5"),
                                            limits = c(.35, 5.15), expand = expansion(mult = 0))
save_figure <- function(plot, name, width, height) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot, width = width, height = height,
         device = cairo_pdf, bg = "white")
  if ("--png" %in% args) ggsave(file.path(out_dir, paste0(name, ".png")), plot,
                               width = width, height = height, dpi = 320, bg = "white")
}

# Mean ARI and the strict acceptable-ARI proportion retain separate y scales/CIs.
for (threshold in c(FALSE, TRUE)) {
  plot <- ggplot(ari_summary, aes(x = delta, color = method, linetype = method, shape = method, group = method))
  if (threshold) plot <- plot + aes(y = prop_gt_070)
  else plot <- plot + aes(y = mean) + geom_hline(yintercept = .70, color = neutral_charcoal,
                                                 linetype = "dashed", linewidth = .55)
  plot <- plot + geom_line(linewidth = .95) + geom_point(size = 2.3)
  if (threshold) plot <- plot + geom_errorbar(aes(ymin = prop_070_ci_low, ymax = prop_070_ci_high), width = .12, linewidth = .65)
  else plot <- plot + geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = .12, linewidth = .65)
  plot <- plot + facet_wrap(~ effect_label, nrow = 1) + method_scales() + delta_scale()
  if (threshold) plot <- plot + scale_y_continuous(breaks = seq(0, 1, .2), labels = percent_format(accuracy = 1), limits = c(0, 1.02))
  else plot <- plot + scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, .2))
  plot <- plot + labs(x = expression(paste("Sampling interval (", Delta, "t)")),
                      y = if (threshold) "Replications above 0.70" else "Adjusted Rand index") + theme_report()
  save_figure(plot, if (threshold) "figure_5_ari_above_070" else "figure_1_ari_means", 8.1, 3.45)
}

heatmap_data <- paired_ari_summary |> mutate(delta_label = factor(
  paste0("Delta t = ", format(delta, nsmall = 1)),
  levels = paste0("Delta t = ", format(c(.5, 1, 5), nsmall = 1))))
plot <- ggplot(heatmap_data, aes(x = delta_label, y = effect_label, fill = mean)) +
  geom_tile(color = "white", linewidth = 1.1) +
  geom_text(aes(label = sprintf("%+.3f", mean)), fontface = "bold",
            color = ifelse(abs(heatmap_data$mean) > .42, "white", "#17202A"), size = 4.1) +
  scale_fill_gradient2(low = aggie_gold, mid = "#F7F7F7", high = aggie_blue, midpoint = 0,
                       limits = c(-.82, .82), labels = number_format(accuracy = .1)) +
  labs(x = NULL, y = "Subgroup-effect magnitude", fill = "CT - DT") + theme_report() +
  theme(legend.position = "right", panel.grid = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text.x = element_text(face = "bold"), axis.text.y = element_text(face = "bold"))
save_figure(plot, "figure_4_ari_difference_heatmap", 7.2, 3.8)

cross_rate_plot_data <- cross_rate_stability_summary |> filter(metric %in% c("stable_correct", "unstable")) |>
  mutate(outcome = factor(recode(as.character(metric), stable_correct = "Correct at every\nsampling interval",
                                 unstable = "Changed across\nsampling intervals"),
                          levels = c("Correct at every\nsampling interval", "Changed across\nsampling intervals")),
         relation_label = factor(recode(as.character(relation),
           `True same-group pairs` = "Pairs from the same\ngenerating subgroup",
           `True different-group pairs` = "Pairs from different\ngenerating subgroups"),
           levels = c("Pairs from the same\ngenerating subgroup", "Pairs from different\ngenerating subgroups")))
plot <- ggplot(cross_rate_plot_data, aes(x = effect, y = mean, color = method, linetype = method, shape = method, group = method)) +
  geom_line(linewidth = .9) + geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = pmax(0, ci_low), ymax = pmin(1, ci_high)), width = .018, linewidth = .62) +
  facet_grid(relation_label ~ outcome) + method_scales(setNames(c("ctS-GIMME", "S-GIMME"), method_levels)) +
  scale_x_continuous(breaks = c(.30, .60, .90), labels = c("0.30", "0.60", "0.90"),
                     limits = c(.25, .95), expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = seq(0, 1, .2), labels = percent_format(accuracy = 1), limits = c(0, 1.01)) +
  labs(x = "Continuous-time subgroup-effect magnitude", y = "Mean proportion of pairs") + theme_report(9.5) +
  theme(strip.text.y = element_text(angle = 0), legend.position = "top")
save_figure(plot, "figure_7_cross_rate_stability", 8.1, 5.15)
message("Wrote the four paper simulation figures to ", normalizePath(out_dir, winslash = "/"))
