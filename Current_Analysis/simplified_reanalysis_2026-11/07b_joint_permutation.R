#!/usr/bin/env Rscript
# T7 part 2. Permutation calibration of the joint tests. Vine labels are permuted
# (8 pseudo-controls among the 40 two-season vines) and dream is refitted with the
# voom weights of the primary fit. Each draw is appended to the output file, so a
# stopped run resumes where it left off.
ROOT_INIT <- Sys.getenv("GRBV_ROOT", if (dir.exists("Current_Analysis")) normalizePath(".") else normalizePath("../.."))
.libPaths(c(file.path(ROOT_INIT, "Reference_Resources/Rlibs_phase2"), .libPaths()))
source(file.path(ROOT_INIT, "Current_Analysis/simplified_reanalysis_2026-11/00_common.R"))
suppressPackageStartupMessages({library(BiocParallel); library(limma); library(variancePartition)})
source(file.path(ROOT, "Current_Analysis/simplified_reanalysis_2026-11/07_shared.R"))
args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args)) as.integer(args[1]) else 200L
vobj <- readRDS(file.path(OUT, "data/vobj_dream.rds"))
ct <- fread(file.path(PRIOR, "Era6_contrast_definitions.csv"))
contrast_list <- setNames(ct$expression, ct$contrast)
form <- ~ 0 + condition + (1 | vine_id)
L <- makeContrastsDream(form, m, contrasts = contrast_list)
two_season <- as.data.table(m)[, .(n = uniqueN(season)), by = vine_id][n == 2]$vine_id
two_season <- sort(as.character(two_season))
stopifnot(length(two_season) == 40)
path <- file.path(OUT, "data/joint_test_permutation_draws.csv")
done <- if (file.exists(path)) fread(path)$draw else integer()
for (b in setdiff(seq_len(B), done)) {
  t0 <- Sys.time()
  set.seed(SEED + 300000L + b)
  pseudo <- sample(two_season, 8)
  m2 <- m
  m2$condition <- make_condition(ifelse(as.character(m$vine_id) %in% pseudo, "neg", "pos"))
  if (nlevels(m2$condition) != 16) stop("Draw ", b, " lost a cell")
  m2$condition <- factor(m2$condition, levels = levels(m$condition))
  fit <- variancePartition::eBayes(dream(vobj, form, m2, L = L, ddf = "Satterthwaite", BPPARAM = BP, quiet = TRUE))
  s <- joint_counts(fit)$summary
  row <- data.table(draw = b, pseudo_controls = paste(sort(pseudo), collapse = ";"),
                    n_real_controls_as_pseudo = sum(pseudo %in% as.character(m$vine_id[m$infection == "neg"])),
                    infection_stage = s$n_fdr05[1], infection_stage_season = s$n_fdr05[2], only_season_specific = s$n_fdr05[3],
                    failures = length(fit$errors), seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  fwrite(row, path, append = file.exists(path))
  cat(sprintf("draw %d: %d / %d / %d (%.0fs)\n", b, row$infection_stage, row$infection_stage_season, row$only_season_specific, row$seconds))
}
d <- fread(path)
obs <- c(infection_stage = 7861, infection_stage_season = 3236, only_season_specific = 1223)
res <- rbindlist(lapply(names(obs), function(n) {
  x <- d[[n]]
  data.table(test = n, observed = obs[[n]], null_median = median(x), null_mean = mean(x), null_p95 = quantile(x, .95), null_max = max(x),
             p95_over_observed = quantile(x, .95) / obs[[n]], perm_p = (1 + sum(x >= obs[[n]])) / (length(x) + 1), n_draws = length(x))
}))
res[, flag_over_10pct := p95_over_observed > .10]
res[, excess_over_p95 := observed - null_p95]
fwrite(res, file.path(OUT, "data/joint_test_permutation_null.csv"))
print(res)
cat("T7 part 2 COMPLETE\n")
