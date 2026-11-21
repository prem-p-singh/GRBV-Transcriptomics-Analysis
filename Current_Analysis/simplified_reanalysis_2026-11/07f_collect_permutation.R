#!/usr/bin/env Rscript
# T7 part 3. Collects the FARM draws and writes the permutation null.
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("GRBV_ROOT", if (dir.exists("Current_Analysis")) normalizePath(".") else normalizePath("../.."))
OUT <- file.path(ROOT, "Current_Analysis/simplified_reanalysis_2026-11")
files <- sort(list.files(file.path(OUT, "data/farm_perm_draws"), pattern = "^draw_\\d+\\.csv$", full.names = TRUE))
d <- rbindlist(lapply(files, fread))[order(draw)]
truth <- d[draw == 0]; d <- d[draw > 0]
# Reconciliation: 200 expected, each draw once.
stopifnot(nrow(truth) == 1, truth$infection_stage == 7861, truth$infection_stage_season == 3236, truth$only_season_specific == 1223)
missing <- setdiff(1:200, d$draw); if (length(missing)) stop("Missing draws: ", paste(missing, collapse = ","))
stopifnot(nrow(d) == 200, !anyDuplicated(d$draw))
fwrite(d, file.path(OUT, "data/joint_test_permutation_draws.csv"))
cat("Draws:", nrow(d), " genes fitted min/median:", min(d$genes_fitted), median(d$genes_fitted),
    " max failures:", max(d$failures), " draws with any failure:", sum(d$failures > 0), "\n")
loc <- fread(file.path(OUT, "data/joint_test_permutation_draws_local_check.csv"))
ck <- merge(loc, d, by = "draw", suffixes = c("_mac", "_farm"))
cat("Mac vs FARM draws compared:", nrow(ck), " same vines:", sum(ck$pseudo_controls_mac == ck$pseudo_controls_farm),
    " same counts:", sum(ck$infection_stage_mac == ck$infection_stage_farm & ck$infection_stage_season_mac == ck$infection_stage_season_farm), "\n")
obs <- c(infection_stage = 7861, infection_stage_season = 3236, only_season_specific = 1223)
res <- rbindlist(lapply(names(obs), function(n) {
  x <- d[[n]]
  data.table(test = n, observed = obs[[n]], null_median = median(x), null_mean = mean(x), null_p95 = unname(quantile(x, .95)),
             null_max = max(x), p95_over_observed = unname(quantile(x, .95)) / obs[[n]], n_draws_nonzero = sum(x > 0),
             perm_p = (1 + sum(x >= obs[[n]])) / (length(x) + 1), n_draws = length(x))
}))
res[, flag_over_10pct := p95_over_observed > .10]
fwrite(res, file.path(OUT, "data/joint_test_permutation_null.csv"))
print(res)
print(d[, .(draws = .N, mean_is = mean(infection_stage), max_is = max(infection_stage), mean_3w = mean(infection_stage_season), max_3w = max(infection_stage_season)), by = n_real_controls_as_pseudo][order(n_real_controls_as_pseudo)])
print(d[order(-infection_stage)][1:5, .(draw, n_real_controls_as_pseudo, infection_stage, infection_stage_season, only_season_specific)])
cat("T7 part 3 COMPLETE\n")
