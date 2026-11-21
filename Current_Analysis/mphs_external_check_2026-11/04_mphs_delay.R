# Step 4: shared delay on the external MPhS scale. Mirrors the in-house displacement():
# per stage x season cell, library-weighted mean of infected minus healthy. The MPhS reference is
# fixed, so draws only resample or relabel vines; nothing is rebuilt.
source("00_common.R")
SEED <- 20261107L; B <- 4999L
SEASONS <- c("2021", "2022")
CELL_KEYS <- as.vector(outer(STAGES, SEASONS, paste, sep = "__"))
cell_of <- factor(paste(audit$stage, audit$season, sep = "__"), CELL_KEYS)
pos <- !audit$healthy

displacement <- function(score, is_pos = pos, wt = rep(1, length(score))) {
  sp <- tapply(score * wt * is_pos, cell_of, sum) / tapply(wt * is_pos, cell_of, sum)
  sn <- tapply(score * wt * !is_pos, cell_of, sum) / tapply(wt * !is_pos, cell_of, sum)
  setNames(as.vector(sp - sn), CELL_KEYS)
}
vine_tab <- unique(audit[, c("vine", "healthy")])
controls <- sort(vine_tab$vine[vine_tab$healthy]); infected <- sort(vine_tab$vine[!vine_tab$healthy])
nseas <- tapply(audit$season, audit$vine, function(z) length(unique(z)))
two_season <- sort(names(nseas)[nseas == 2])
stopifnot(length(controls) == 8, length(infected) == 36, length(two_season) == 40)

infer <- function(score, label) {
  est <- displacement(score)
  boot <- t(vapply(seq_len(B), function(b) {
    set.seed(SEED + b)
    m <- c(table(sample(controls, 8, TRUE)), table(sample(infected, 36, TRUE)))
    wt <- unname(m[audit$vine]); wt[is.na(wt)] <- 0
    displacement(score, wt = wt)
  }, numeric(8)))
  perm <- t(vapply(seq_len(B), function(b) {
    set.seed(SEED + 100000L + b)
    pseudo <- sample(two_season, 8)
    displacement(score, is_pos = !(audit$vine %in% pseudo))
  }, numeric(8)))
  res <- data.frame(score = label, cell = CELL_KEYS, estimate = est,
    ci_low = apply(boot, 2, quantile, .025, na.rm = TRUE), ci_high = apply(boot, 2, quantile, .975, na.rm = TRUE),
    n_boot = colSums(is.finite(boot)),
    perm_p = (1 + colSums(abs(perm) >= rep(abs(est), each = B))) / (B + 1), n_perm = B)
  res$perm_fdr <- p.adjust(res$perm_p, "BH")
  d <- boot[, "Harvest__2021"] - boot[, "Harvest__2022"]
  sd <- data.frame(score = label, contrast = "Harvest 2021 minus Harvest 2022",
    estimate = est[["Harvest__2021"]] - est[["Harvest__2022"]],
    ci_low = quantile(d, .025, na.rm = TRUE), ci_high = quantile(d, .975, na.rm = TRUE), n_boot = sum(is.finite(d)))
  list(cells = res, season = sd)
}

x_cpm <- expr("cpm")
audit$mphs <- project(x_cpm, audit$Sample, "scale")
ind_genes <- gmap[gmap$inhouse_marker != "TRUE", ]
stopifnot(nrow(ind_genes) >= 120)
audit$mphs_indep <- project(expr("cpm", ind_genes), audit$Sample, "scale")

inh <- read.csv("../simplified_reanalysis_2026-11/data/dev_score_predictions.csv", colClasses = c(Sample = "character"))
audit$new_score <- inh$new_score[match(audit$Sample, inh$Sample)]
stopifnot(!anyNA(audit$new_score))
write.csv(audit[, c("Sample", "vine", "season", "stage", "stage_num", "infection", "mphs", "mphs_indep", "new_score")],
          "data/mphs_all_projection.csv", row.names = FALSE)

main <- infer(audit$mphs, "mphs_primary"); ind <- infer(audit$mphs_indep, "mphs_non_marker_genes")
write.csv(main$cells, "data/mphs_delay_inference.csv", row.names = FALSE)
write.csv(rbind(main$season, ind$season), "data/mphs_season_difference.csv", row.names = FALSE)
write.csv(ind$cells, "data/mphs_independent_genes_delay.csv", row.names = FALSE)

# Agreement with the in-house score
sp <- function(i) cor(audit$mphs[i], audit$new_score[i], method = "spearman")
fit <- lm(new_score ~ mphs, data = audit[audit$healthy, ])
conv <- coef(fit)[["mphs"]]
inh_cells <- read.csv("../simplified_reanalysis_2026-11/data/dev_displacement_inference.csv")
inh_cells <- inh_cells[match(CELL_KEYS, inh_cells$cell), ]
agree <- data.frame(cell = CELL_KEYS, mphs_delay = main$cells$estimate,
  mphs_delay_on_inhouse_scale = main$cells$estimate * conv, inhouse_delay = inh_cells$estimate)
agree$sign_agree <- sign(agree$mphs_delay) == sign(agree$inhouse_delay)
agree$ratio_converted_over_inhouse <- agree$mphs_delay_on_inhouse_scale / agree$inhouse_delay
lib <- data.frame(cell = c("library_rho_all327", "library_rho_healthy", "library_rho_infected",
                           "healthy_fit_slope_inhouse_per_mphs_step", "healthy_fit_r2",
                           "harvest_2021_delay_larger_than_2022_mphs", "harvest_2021_delay_larger_than_2022_mphs_non_marker"),
  mphs_delay = c(sp(TRUE), sp(audit$healthy), sp(!audit$healthy), conv, summary(fit)$r.squared,
                 main$cells$estimate[4] < main$cells$estimate[8], ind$cells$estimate[4] < ind$cells$estimate[8]))
agree_out <- merge(agree, lib, all = TRUE, sort = FALSE)
write.csv(agree_out, "data/mphs_vs_inhouse_agreement.csv", row.names = FALSE)

# Figure: (A) per-cell delay with CIs, (B) per-library MPhS vs in-house score
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
cmp <- rbind(
  data.frame(scale = "MPhS (converted to in-house units)", cell = CELL_KEYS, est = main$cells$estimate * conv,
             lo = main$cells$ci_low * conv, hi = main$cells$ci_high * conv),
  data.frame(scale = "In-house healthy-only score", cell = CELL_KEYS, est = inh_cells$estimate,
             lo = inh_cells$nested_ci_low, hi = inh_cells$nested_ci_high))
cmp$stage <- factor(sub("__.*", "", cmp$cell), STAGES); cmp$season <- sub(".*__", "", cmp$cell)
pa <- ggplot(cmp, aes(stage, est, colour = scale)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.5), size = 0.35) +
  facet_wrap(~season) + scale_colour_manual(values = c("#7F2944", "#087F8C"), name = NULL) +
  labs(x = NULL, y = "Infected minus healthy (in-house stage units)", title = "A. Shared delay, 95% vine-bootstrap CIs") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))
pb <- ggplot(audit, aes(new_score, mphs, colour = ifelse(healthy, "Healthy", "Infected"))) +
  geom_point(position = position_jitter(height = 0.15, width = 0, seed = 1), size = 1, alpha = 0.6) +
  scale_colour_manual(values = c(Healthy = "#087F8C", Infected = "#7F2944"), name = NULL) +
  labs(x = "In-house developmental score", y = "MPhS stage (1-30)",
       title = sprintf("B. Per library (Spearman rho = %.2f)", sp(TRUE))) + theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave("figures/mphs_delay_comparison.png", pa + pb + plot_layout(widths = c(1.6, 1)), width = 11, height = 4.8, dpi = 300)

writeLines(sprintf("seed %d; bootstrap and permutation draws %d each; draw b uses seed+b (bootstrap) and seed+100000+b (permutation)", SEED, B), "logs/04_seeds.txt")
print(main$cells, digits = 3); print(ind$cells[, c("cell", "estimate", "ci_low", "ci_high", "perm_fdr")], digits = 3)
print(rbind(main$season, ind$season), digits = 3); print(agree_out, digits = 3)
