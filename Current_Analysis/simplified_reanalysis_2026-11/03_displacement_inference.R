#!/usr/bin/env Rscript
# T3. Nested vine bootstrap and vine-label permutation for the shared delay.
# Both rerun the full chain (markers, centroids, projection) inside each draw,
# with leave-one-reference-vine-out cross-fitting, so reference uncertainty is
# included. Each draw has its own seed, so results do not depend on core count.
ROOT_INIT <- Sys.getenv("GRBV_ROOT", if (dir.exists("Current_Analysis")) normalizePath(".") else normalizePath("../.."))
.libPaths(c(file.path(ROOT_INIT, "Reference_Resources/Rlibs_phase2"), .libPaths()))
source(file.path(ROOT_INIT, "Current_Analysis/simplified_reanalysis_2026-11/00_common.R"))
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
B_BOOT <- if (length(args) >= 1) as.integer(args[1]) else 1999L
B_PERM <- if (length(args) >= 2) as.integer(args[2]) else 4999L
smoke <- B_BOOT < 1999L || B_PERM < 4999L

vines <- unique(meta[, .(vine_id, infection)])
controls <- sort(vines[infection == "neg"]$vine_id)
infected <- sort(vines[infection == "pos"]$vine_id)
seasons_per_vine <- meta[, .(n = uniqueN(season)), by = vine_id]
two_season <- sort(seasons_per_vine[n == 2]$vine_id)
one_season <- setdiff(meta$vine_id, two_season)
stopifnot(length(controls) == 8, length(infected) == 36, length(two_season) == 40,
          length(one_season) == 4, all(one_season %in% infected))

# Cross-fitted score. ref_vines: unique reference vines; w: their multiplicities;
# other: remaining vines, spread at random over the folds.
crossfit <- function(ref_vines, w, other) {
  m <- length(ref_vines)
  fold_of <- setNames(c(seq_len(m), sample(rep_len(sample(m), length(other)))), c(ref_vines, other))
  score <- rep(NA_real_, ncol(E))
  for (f in seq_len(m)) {
    wf <- w; wf[f] <- 0
    ref <- make_ref(setNames(wf, ref_vines))
    ids <- pick_markers(ref, background)
    idx <- which(fold_of[meta$vine_id] == f)
    score[idx] <- predict_position(ref, ids, idx)
  }
  score
}

boot_draw <- function(b) {
  set.seed(SEED + b)
  ctr <- table(sample(controls, length(controls), replace = TRUE))
  inf <- table(sample(infected, length(infected), replace = TRUE))
  if (length(ctr) < 2) return(rep(NA_real_, 8))  # cannot cross-fit one vine
  mult <- setNames(numeric(nrow(vines)), vines$vine_id)
  mult[names(ctr)] <- as.numeric(ctr); mult[names(inf)] <- as.numeric(inf)
  score <- crossfit(names(ctr), as.numeric(ctr), names(inf))
  wt <- unname(mult[meta$vine_id])
  keep <- wt > 0
  score[!keep] <- 0
  suppressWarnings(displacement(score, wt = wt))
}

perm_draw <- function(b) {
  set.seed(SEED + 100000L + b)
  pseudo <- sort(sample(two_season, length(controls)))
  labels <- setNames(rep("pos", nrow(vines)), vines$vine_id)
  labels[pseudo] <- "neg"
  score <- crossfit(pseudo, rep(1, length(pseudo)), setdiff(vines$vine_id, pseudo))
  displacement(score, labels)
}

t0 <- Sys.time()
boot <- do.call(rbind, mclapply(seq_len(B_BOOT), boot_draw, mc.cores = CORES, mc.preschedule = TRUE))
cat("Bootstrap draws:", nrow(boot), "seconds:", round(difftime(Sys.time(), t0, units = "secs")), "\n")
t1 <- Sys.time()
perm <- do.call(rbind, mclapply(seq_len(B_PERM), perm_draw, mc.cores = CORES, mc.preschedule = TRUE))
cat("Permutation draws:", nrow(perm), "seconds:", round(difftime(Sys.time(), t1, units = "secs")), "\n")
stopifnot(is.numeric(boot), is.numeric(perm), ncol(boot) == 8, ncol(perm) == 8, !anyNA(perm))

pred <- fread(file.path(OUT, "data/dev_score_predictions.csv"))
pred <- pred[match(meta$Sample, as.character(Sample))]
est <- displacement(pred$new_score)
old <- fread(file.path(PRIOR, "Era7_global_stage_displacement.csv"))
old <- old[match(CELL_KEYS, cell)]

res <- data.table(cell = CELL_KEYS, estimate = est, old_estimate = old$effect,
                  old_ci_low = old$ci_low, old_ci_high = old$ci_high, old_fdr = old$FDR,
                  nested_ci_low = apply(boot, 2, quantile, .025, na.rm = TRUE),
                  nested_ci_high = apply(boot, 2, quantile, .975, na.rm = TRUE),
                  n_boot_used = colSums(is.finite(boot)),
                  perm_p = (1 + colSums(abs(perm) >= rep(abs(est), each = nrow(perm)))) / (nrow(perm) + 1),
                  perm_null_sd = apply(perm, 2, sd), n_perm = nrow(perm))
res[, perm_fdr := p.adjust(perm_p, "BH")]
res[, c("stage", "season") := tstrsplit(cell, "__", fixed = TRUE)]

# Season difference at harvest: is the 2021 delay larger than the 2022 delay?
h21 <- "Harvest__2021"; h22 <- "Harvest__2022"
d_est <- est[h21] - est[h22]
d_boot <- boot[, h21] - boot[, h22]
d_perm <- perm[, h21] - perm[, h22]
season_diff <- data.table(contrast = "Harvest 2021 minus Harvest 2022", estimate = d_est,
                          nested_ci_low = quantile(d_boot, .025, na.rm = TRUE),
                          nested_ci_high = quantile(d_boot, .975, na.rm = TRUE),
                          perm_p = (1 + sum(abs(d_perm) >= abs(d_est))) / (length(d_perm) + 1),
                          n_boot_used = sum(is.finite(d_boot)), n_perm = length(d_perm))
suffix <- if (smoke) "_SMOKE" else ""
fwrite(res, file.path(OUT, sprintf("data/dev_displacement_inference%s.csv", suffix)))
fwrite(season_diff, file.path(OUT, sprintf("data/dev_displacement_season_difference%s.csv", suffix)))
if (!smoke) {
  fwrite(as.data.table(boot), file.path(OUT, "data/dev_displacement_bootstrap_draws.csv.gz"))
  fwrite(as.data.table(perm), file.path(OUT, "data/dev_displacement_permutation_draws.csv.gz"))
}
print(res[, .(cell, estimate, nested_ci_low, nested_ci_high, perm_p, perm_fdr, n_boot_used)])
print(season_diff)
cat("T3 COMPLETE\n")
