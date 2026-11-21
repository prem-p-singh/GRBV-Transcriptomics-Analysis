#!/usr/bin/env Rscript
# T8. Within-infected gradients (prespecified secondary, nulls included). Uses viral
# reads from the same library and infection cohort. The qPCR titer is not used.
ROOT_INIT <- Sys.getenv("GRBV_ROOT", if (dir.exists("Current_Analysis")) normalizePath(".") else normalizePath("../.."))
.libPaths(c(file.path(ROOT_INIT, "Reference_Resources/Rlibs_phase2"), .libPaths()))
source(file.path(ROOT_INIT, "Current_Analysis/simplified_reanalysis_2026-11/00_common.R"))
suppressPackageStartupMessages({library(lmerTest)})
audit <- fread(file.path(OUT, "data/sample_audit.csv"))[match(meta$Sample, as.character(Sample))]
pred <- fread(file.path(OUT, "data/dev_score_predictions.csv"))[match(meta$Sample, as.character(Sample))]
R <- readRDS(file.path(OUT, "data/dream_fits.rds"))$residuals
stopifnot(identical(colnames(R), meta$Sample))

controls <- sort(unique(meta[infection == "neg"]$vine_id))
ps <- intersect(mem[set_id == "BRANCH:58"]$gene, rownames(E))
branch <- rep(NA_real_, ncol(E))
for (f in sort(unique(meta$fold))) {
  drop <- unique(meta[infection == "neg" & fold == f]$vine_id)
  ref <- make_ref(setNames(as.numeric(!controls %in% drop), controls))
  idx <- which(meta$fold == f)
  branch[idx] <- predict_position(ref, pick_markers(ref, ps, 100L), idx)
}
shsp <- intersect(mem[set_id == "BRANCH:50"]$gene, rownames(E))
d <- data.table(Sample = meta$Sample, vine = meta$vine_id, infection = meta$infection, stage = meta$Maturity_stage, season = meta$season,
                cohort = audit$cohort, log10_viral = log10(audit$viral_cpm),
                developmental_offset = pred$new_score - pred$stage_number,
                photosystem_offset = branch - pred$new_score,
                shsp_adjusted_residual = colMeans(R[shsp, , drop = FALSE]))[infection == "pos"]
stopifnot(nrow(d) == 264, !anyNA(d))
fwrite(d, file.path(OUT, "data/within_infected_library_values.csv"))
d[, cell := interaction(stage, season, drop = TRUE)]
print(d[, .(libraries = .N, vines = uniqueN(vine)), by = cohort][order(cohort)])
outcomes <- c("developmental_offset", "photosystem_offset", "shsp_adjusted_residual")
res <- list()
for (y in outcomes) {
  f <- lmer(as.formula(paste(y, "~ log10_viral + cell + (1 | vine)")), d, REML = TRUE)
  s <- coef(summary(f))["log10_viral", ]
  res[[length(res) + 1]] <- data.table(outcome = y, predictor = "log10 viral CPM", term = "slope per log10", estimate = s["Estimate"],
    se = s["Std. Error"], df = s["df"], statistic = s["t value"], p = s["Pr(>|t|)"], n_libraries = nrow(d), n_vines = uniqueN(d$vine), singular = isSingular(f))
  dc <- d[cohort %in% c("2016", "2019", "2020", "2021")]
  dc[, cohort_f := factor(cohort, levels = c("2016", "2019", "2020", "2021"))]
  dc[, years := as.numeric(cohort)]
  g <- lmer(as.formula(paste(y, "~ cohort_f + cell + (1 | vine)")), dc, REML = TRUE)
  a <- anova(g)["cohort_f", ]
  res[[length(res) + 1]] <- data.table(outcome = y, predictor = "infection cohort", term = "joint F, 3 df", estimate = NA_real_,
    se = NA_real_, df = a[["DenDF"]], statistic = a[["F value"]], p = a[["Pr(>F)"]], n_libraries = nrow(dc), n_vines = uniqueN(dc$vine), singular = isSingular(g))
  h <- lmer(as.formula(paste(y, "~ years + cell + (1 | vine)")), dc, REML = TRUE)
  s <- coef(summary(h))["years", ]
  res[[length(res) + 1]] <- data.table(outcome = y, predictor = "infection cohort", term = "linear trend per cohort year", estimate = s["Estimate"],
    se = s["Std. Error"], df = s["df"], statistic = s["t value"], p = s["Pr(>|t|)"], n_libraries = nrow(dc), n_vines = uniqueN(dc$vine), singular = isSingular(h))
}
res <- rbindlist(res)
res[, fdr := p.adjust(p, "BH")]
fwrite(res, file.path(OUT, "data/within_infected_gradients.csv"))
print(res)
cat("T8 COMPLETE\n")
