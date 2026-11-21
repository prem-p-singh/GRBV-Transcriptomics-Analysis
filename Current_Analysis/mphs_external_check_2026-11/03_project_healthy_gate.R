# Step 3: project the 63 healthy libraries and apply the pre-registered gate (logs/03_gate_thresholds.md).
source("00_common.R"); set.seed(20261107)
h <- audit$Sample[audit$healthy]
x_cpm <- expr("cpm"); x_tmm <- expr("tmm_logcpm")
p <- data.frame(Sample = h,
  mphs_primary = project(x_cpm, h, "scale"),
  mphs_alt1_tmm_logcpm = project(x_tmm, h, "scale"),
  mphs_alt2_cpm_nolength = project(x_cpm, h, "scale"),
  mphs_extra_means_SDs = project(x_cpm, h, "means_SDs"))
all_primary <- project(x_cpm, audit$Sample, "scale")      # single all-327 projection (step 4 scaling)
p$mphs_from_all327 <- all_primary[h]
p <- merge(audit[audit$healthy, c("Sample", "vine", "season", "stage", "stage_num")], p, by = "Sample")
write.csv(p, "data/mphs_healthy_projection.csv", row.names = FALSE)

gate_stats <- function(v, d = p) {
  do.call(rbind, lapply(split(d, d$season), function(s) {
    mu <- tapply(s[[v]], s$stage_num, mean); sdw <- sqrt(sum(tapply(s[[v]], s$stage_num, function(z) sum((z - mean(z))^2))) / (nrow(s) - 4))
    data.frame(score = v, season = s$season[1], n = nrow(s),
      mean_pre = mu[1], mean_ver = mu[2], mean_post = mu[3], mean_harv = mu[4],
      monotonic = all(diff(mu) > 0),
      spearman_rho = suppressWarnings(cor(s[[v]], s$stage_num, method = "spearman")),
      pooled_within_sd = sdw, n_pairs_separated = sum(diff(mu) > sdw),
      harvest_frac_within1_of_max = mean(s[[v]][s$stage_num == 3] >= 29))
  }))
}
res <- do.call(rbind, lapply(c("mphs_primary", "mphs_alt1_tmm_logcpm", "mphs_extra_means_SDs", "mphs_from_all327"), gate_stats))
rob <- c(alt1 = cor(p$mphs_primary, p$mphs_alt1_tmm_logcpm, method = "spearman"),
         alt2 = cor(p$mphs_primary, p$mphs_alt2_cpm_nolength, method = "spearman"),
         extra_means_SDs = cor(p$mphs_primary, p$mphs_extra_means_SDs, method = "spearman"),
         from_all327 = cor(p$mphs_primary, p$mphs_from_all327, method = "spearman"))
pr <- res[res$score == "mphs_primary", ]
crit <- data.frame(
  criterion = c("monotonic_both_seasons", "spearman_ge_0.85_each_season", "separation_2of3_each_season",
                "robust_alt1_ge_0.90", "robust_alt2_ge_0.90"),
  value = c(paste(pr$monotonic, collapse = "/"), paste(round(pr$spearman_rho, 3), collapse = "/"),
            paste(pr$n_pairs_separated, collapse = "/"), round(rob["alt1"], 3), round(rob["alt2"], 3)),
  pass = c(all(pr$monotonic), all(pr$spearman_rho >= 0.85), all(pr$n_pairs_separated >= 2),
           rob["alt1"] >= 0.90, rob["alt2"] >= 0.90))
crit <- rbind(crit, data.frame(criterion = "GATE_PASSED", value = "", pass = all(crit$pass)))
write.csv(res, "data/mphs_gate_stage_summary.csv", row.names = FALSE)
write.csv(rbind(crit, data.frame(criterion = paste0("report_only_rho_primary_vs_", names(rob)[3:4]), value = round(rob[3:4], 3), pass = NA)),
          "data/mphs_gate_results.csv", row.names = FALSE)
rng <- aggregate(mphs_primary ~ season + stage_num, p, function(z) paste(min(z), max(z), sep = "-"))
write.csv(rng, "data/mphs_healthy_stage_ranges.csv", row.names = FALSE)

suppressPackageStartupMessages(library(ggplot2))
p$stage <- factor(p$stage, STAGES)
g <- ggplot(p, aes(stage, mphs_primary, colour = factor(season))) +
  geom_line(aes(group = interaction(vine, season)), alpha = 0.35) +
  geom_point(position = position_jitter(width = 0.08, height = 0.15, seed = 1), size = 1.6) +
  stat_summary(aes(group = season), fun = mean, geom = "line", linewidth = 1.1) +
  scale_colour_manual(values = c(`2021` = "#087F8C", `2022` = "#7F2944"), name = "Season") +
  labs(x = "Nominal stage", y = "MPhS stage (1-30)", title = "Healthy libraries on the MPhS scale (primary normalization)") +
  theme_bw(base_size = 11)
ggsave("figures/mphs_healthy_gate.png", g, width = 7, height = 4.5, dpi = 300)
print(res, digits = 3); print(crit); print(rob); print(rng)
if (!crit$pass[crit$criterion == "GATE_PASSED"]) writeLines(c("# Gate failed", "", capture.output(print(crit))), "logs/gate_failed.md")
