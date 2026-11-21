#!/usr/bin/env Rscript
# T2. Uncensored developmental score on one pooled-season scale.
ROOT_INIT <- Sys.getenv("GRBV_ROOT", if (dir.exists("Current_Analysis")) normalizePath(".") else normalizePath("../.."))
.libPaths(c(file.path(ROOT_INIT, "Reference_Resources/Rlibs_phase2"), .libPaths()))
source(file.path(ROOT_INIT, "Current_Analysis/simplified_reanalysis_2026-11/00_common.R"))
set.seed(SEED)
controls <- sort(unique(meta[infection == "neg"]$vine_id))
stopifnot(length(controls) == 8)
fold_vine <- unique(meta[infection == "neg", .(vine_id, fold)])
weights_without <- function(drop) setNames(as.numeric(!controls %in% drop), controls)

# Smooth sensitivity: natural cubic spline per marker through the pooled centroids.
GRID <- seq(-1.5, 4.5, by = 0.01)
predict_spline <- function(ref, ids, idx) {
  curve <- t(apply(ref$pooled[ids, , drop = FALSE], 1, function(y) splinefun(0:3, y, method = "natural")(GRID))) / ref$scale[ids]
  obs <- E[ids, idx, drop = FALSE] / ref$scale[ids]
  d <- outer(colSums(curve^2), colSums(obs^2), "+") - 2 * crossprod(curve, obs)
  GRID[apply(d, 2, which.min)]
}

variants <- list(new = list(pool = TRUE, extrapolate = TRUE),
                 seasonal_extrap = list(pool = FALSE, extrapolate = TRUE),
                 pooled_clamped = list(pool = TRUE, extrapolate = FALSE),
                 repro_old = list(pool = FALSE, extrapolate = FALSE))
scores <- sapply(names(variants), function(n) rep(NA_real_, ncol(E)))
spline <- rep(NA_real_, ncol(E))
for (f in fold_vine$fold) {
  drop <- fold_vine[fold == f]$vine_id
  idx <- which(meta$fold == f)
  for (n in names(variants)) {
    v <- variants[[n]]
    ref <- make_ref(weights_without(drop), pool = v$pool)
    scores[idx, n] <- predict_position(ref, pick_markers(ref, background), idx, v$extrapolate)
  }
  ref <- make_ref(weights_without(drop), pool = TRUE)
  spline[idx] <- predict_spline(ref, pick_markers(ref, background), idx)
}
stopifnot(!anyNA(scores), !anyNA(spline))

# Regression guard: the old settings must reproduce the Era 7 scores exactly.
old <- fread(file.path(PRIOR, "Era7_developmental_predictions.csv"))[model == "background_1000"]
old <- old[match(meta$Sample, as.character(Sample))]
stopifnot(!anyNA(old$predicted_stage))
worst <- max(abs(old$predicted_stage - scores[, "repro_old"]))
cat("Max |repro - Era 7 score|:", worst, "\n")
if (worst > 1e-8) stop("Refactored scoring does not reproduce Era 7")

pred <- meta[, .(Sample, vine_id, season, stage = Maturity_stage, stage_number, infection, fold)]
pred[, `:=`(old_score = old$predicted_stage, new_score = scores[, "new"],
            seasonal_extrap = scores[, "seasonal_extrap"], pooled_clamped = scores[, "pooled_clamped"],
            spline_score = spline)]
fwrite(pred, file.path(OUT, "data/dev_score_predictions.csv"))

# Held-out season: the reference never sees the target season.
hold <- rbindlist(lapply(SEASONS, function(y) {
  idx <- which(meta$infection == "neg" & meta$season == y)
  train_ref <- make_ref(setNames(rep(1, length(controls)), controls), pool = TRUE,
                        allowed = meta$season != y)
  ids <- pick_markers(train_ref, background)
  old_s <- fread(file.path(PRIOR, "Era7_developmental_season_holdout.csv"))[model == "background_1000"]
  z <- meta[idx, .(Sample, vine_id, season, stage_number)]
  z[, `:=`(new_score = predict_position(train_ref, ids, idx),
           old_score = old_s$predicted_stage[match(Sample, as.character(old_s$Sample))])]
  z
}))
stopifnot(!anyNA(hold$old_score))
fwrite(hold, file.path(OUT, "data/dev_score_heldout_season.csv"))
nearest <- function(x) pmin(3, pmax(0, round(x)))
val <- function(d, label) rbindlist(lapply(c("old_score", "new_score"), function(sc) {
  x <- d[[sc]]
  data.table(validation = label, score = sc, n = nrow(d), RMSE = sqrt(mean((x - d$stage_number)^2)),
             MAE = mean(abs(x - d$stage_number)), spearman = cor(x, d$stage_number, method = "spearman"),
             nearest_correct = sum(nearest(x) == d$stage_number), nearest_accuracy = mean(nearest(x) == d$stage_number))
}))
valid <- rbind(val(pred[infection == "neg"], "held_out_vine"), val(hold, "held_out_season"))
fwrite(valid, file.path(OUT, "data/dev_score_validation.csv"))
print(valid)

long <- melt(pred, id.vars = c("infection", "season", "stage", "stage_number"),
             measure.vars = c("old_score", "new_score"), variable.name = "score", value.name = "x")
bounds <- long[, .(n = .N, at_or_below_0 = sum(x <= 0), at_or_above_3 = sum(x >= 3),
                   exactly_0 = sum(x == 0), exactly_3 = sum(x == 3),
                   frac_at_or_below_0 = mean(x <= 0), frac_at_or_above_3 = mean(x >= 3),
                   min_score = min(x), max_score = max(x)),
               by = .(score, infection, season, stage)]
fwrite(bounds[order(score, infection, season, match(stage, STAGES))], file.path(OUT, "data/dev_score_bounds.csv"))
print(bounds[score == "old_score" & (exactly_0 > 0 | exactly_3 > 0)])
print(bounds[score == "new_score", .(exactly_0 = sum(exactly_0), exactly_3 = sum(exactly_3), min = min(min_score), max = max(max_score))])

disp <- data.table(cell = CELL_KEYS)
for (cl in c("old_score", "new_score", "seasonal_extrap", "pooled_clamped", "spline_score"))
  disp[[cl]] <- displacement(pred[[cl]])
disp[, c("stage", "season") := tstrsplit(cell, "__", fixed = TRUE)]
ctr <- pred[infection == "neg", .(control_mean_new = mean(new_score), n_control = .N), by = .(stage, season)]
inf <- pred[infection == "pos", .(infected_mean_new = mean(new_score), n_infected = .N), by = .(stage, season)]
disp <- merge(merge(disp, ctr, by = c("stage", "season")), inf, by = c("stage", "season"))
disp <- disp[order(season, match(stage, STAGES))]
fwrite(disp, file.path(OUT, "data/dev_displacement.csv"))
print(disp[, .(cell, old_score, new_score, seasonal_extrap, pooled_clamped, spline_score)])
h <- disp[stage == "Harvest"]
cat(sprintf("Harvest displacement 2021 %.3f, 2022 %.3f; 2021 larger in magnitude: %s\n",
            h[season == "2021"]$new_score, h[season == "2022"]$new_score,
            abs(h[season == "2021"]$new_score) > abs(h[season == "2022"]$new_score)))
cat("T2 COMPLETE\n")
