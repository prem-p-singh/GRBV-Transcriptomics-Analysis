#!/usr/bin/env Rscript
# T4. Development-adjusted expression. The adjusted effect is the infection
# effect not accounted for by the shared developmental position. It is not a
# causal direct effect, because development is itself affected by infection.
ROOT_INIT <- Sys.getenv("GRBV_ROOT", if (dir.exists("Current_Analysis")) normalizePath(".") else normalizePath("../.."))
.libPaths(c(file.path(ROOT_INIT, "Reference_Resources/Rlibs_phase2"), .libPaths()))
source(file.path(ROOT_INIT, "Current_Analysis/simplified_reanalysis_2026-11/00_common.R"))
suppressPackageStartupMessages({library(BiocParallel); library(edgeR); library(limma); library(variancePartition)})
STAGE_CODES <- c(Pre_veraison = "Pre", Veraison = "Veraison", Post_veraison = "Post", Harvest = "Harvest")
FDR <- 0.05
BP <- MulticoreParam(CORES, progressbar = FALSE)

m <- as.data.frame(meta)
rownames(m) <- m$Sample
m$condition <- factor(paste(m$infection, STAGE_CODES[as.character(m$Maturity_stage)], m$season, sep = "_"))
m$vine_id <- factor(m$vine_id)
m$stage <- factor(as.character(m$Maturity_stage), levels = STAGES)
counts <- obj$counts
stopifnot(identical(colnames(counts), m$Sample), identical(rownames(counts), rownames(E)))
form <- ~ 0 + condition + (1 | vine_id)

cell_expr <- function(s, y) sprintf("conditionpos_%s_%s - conditionneg_%s_%s", STAGE_CODES[[s]], y, STAGE_CODES[[s]], y)
contrasts_all <- c()
for (y in SEASONS) for (s in STAGES) contrasts_all[paste("infection", STAGE_CODES[[s]], y, sep = "__")] <- cell_expr(s, y)
L <- makeContrastsDream(form, m, contrasts = contrasts_all)

cache <- file.path(OUT, "data/vobj_dream.rds")
if (file.exists(cache)) vobj <- readRDS(cache) else {
  d <- calcNormFactors(DGEList(counts = counts), method = "TMM")
  vobj <- voomWithDreamWeights(d, form, m, plot = FALSE, BPPARAM = BP)
  saveRDS(vobj, cache)
}
cat("voom E vs Era 7 E, max abs difference:", max(abs(vobj$E - E)), "\n")
stopifnot(max(abs(vobj$E - E)) < 1e-8)

extract <- function(fit, label) rbindlist(lapply(names(contrasts_all), function(cn) {
  tt <- variancePartition::topTable(fit, coef = cn, number = Inf, sort.by = "none", adjust.method = "BH")
  f <- strsplit(cn, "__", fixed = TRUE)[[1]]
  data.table(gene = rownames(tt), model = label, stage = names(STAGE_CODES)[match(f[2], STAGE_CODES)], season = f[3],
             log2FC = tt$logFC, t = tt$t, z_std = tt$z.std, p_value = tt$P.Value, fdr = tt$adj.P.Val)
}))
fit_model <- function(v, f, L, label) {
  t0 <- Sys.time()
  fit <- dream(v, f, m, L = L, ddf = "Satterthwaite", BPPARAM = BP, quiet = TRUE)
  fit <- variancePartition::eBayes(fit)
  cat(label, "seconds:", round(difftime(Sys.time(), t0, units = "secs")), "failures:", length(fit$errors), "\n")
  list(fit = fit, table = extract(fit, label))
}

# Unadjusted refit: must match the primary Era 6 results.
primary <- fit_model(vobj, form, L, "unadjusted")
era6 <- fread(file.path(PRIOR, "Era6_cell_contrasts.csv"))
chk <- merge(primary$table, era6[, .(gene, stage, season = as.character(season), e6_lfc = log2FC, e6_fdr = fdr_within_contrast)],
             by = c("gene", "stage", "season"))
stopifnot(nrow(chk) == 14604 * 8)
cat("Max |log2FC - Era 6|:", max(abs(chk$log2FC - chk$e6_lfc)),
    " FDR<0.05 mine:", sum(chk$fdr < FDR), " Era 6:", sum(chk$e6_fdr < FDR), "\n")
stopifnot(cor(chk$log2FC, chk$e6_lfc) > 0.99999)

# Primary adjustment: subtract expression expected at each library's own score,
# on the cross-fitted reference from that library's fold.
pred <- fread(file.path(OUT, "data/dev_score_predictions.csv"))
pred <- pred[match(meta$Sample, as.character(Sample))]
stopifnot(!anyNA(pred$new_score))
controls <- sort(unique(meta[infection == "neg"]$vine_id))
R <- E * NA_real_
for (f in sort(unique(meta$fold))) {
  drop <- unique(meta[infection == "neg" & fold == f]$vine_id)
  ref <- make_ref(setNames(as.numeric(!controls %in% drop), controls))
  idx <- which(meta$fold == f)
  # Guard: the stored score must be the one this reference gives.
  stopifnot(max(abs(predict_position(ref, pick_markers(ref, background), idx) - pred$new_score[idx])) < 1e-8)
  R[, idx] <- E[, idx, drop = FALSE] - expected_expression(ref, rownames(E), pred$new_score[idx])
}
stopifnot(!anyNA(R))
vres <- vobj; vres$E <- R
adjusted <- fit_model(vres, form, L, "adjusted_residual")

# Sensitivity: score offset as a covariate with stage-specific slopes.
m$offset_score <- pred$new_score - pred$stage_number
form_cov <- ~ 0 + condition + offset_score:stage + (1 | vine_id)
L_cov <- makeContrastsDream(form_cov, m, contrasts = contrasts_all)
covariate <- fit_model(vobj, form_cov, L_cov, "adjusted_covariate")

saveRDS(list(primary = primary$fit, adjusted = adjusted$fit, residuals = R), file.path(OUT, "data/dream_fits.rds"))

tabs <- rbindlist(list(primary$table, adjusted$table, covariate$table))
fwrite(tabs[model != "unadjusted"], file.path(OUT, "data/adjusted_dream_cell_contrasts.csv"))
fwrite(primary$table, file.path(OUT, "data/unadjusted_dream_cell_contrasts.csv"))

w <- dcast(tabs, gene + stage + season ~ model, value.var = c("log2FC", "t", "fdr"))
agree <- w[, .(spearman_t = cor(t_adjusted_residual, t_adjusted_covariate, method = "spearman"),
               pearson_log2FC = cor(log2FC_adjusted_residual, log2FC_adjusted_covariate),
               n_residual = sum(fdr_adjusted_residual < FDR), n_covariate = sum(fdr_adjusted_covariate < FDR),
               n_both = sum(fdr_adjusted_residual < FDR & fdr_adjusted_covariate < FDR)), by = .(season, stage)]
agree[, jaccard := n_both / (n_residual + n_covariate - n_both)]
fwrite(agree, file.path(OUT, "data/adjustment_method_agreement.csv"))

w[, class := fifelse(fdr_adjusted_residual < FDR, "beyond_delay",
                     fifelse(fdr_unadjusted < FDR, "explained_by_delay", "not_changed"))]
w[, `:=`(revealed_by_adjustment = fdr_adjusted_residual < FDR & fdr_unadjusted >= FDR,
         sign_flip = fdr_adjusted_residual < FDR & fdr_unadjusted < FDR &
           sign(log2FC_adjusted_residual) != sign(log2FC_unadjusted))]
markers_any <- rbindlist(lapply(sort(unique(meta$fold)), function(f) {
  drop <- unique(meta[infection == "neg" & fold == f]$vine_id)
  data.table(gene = pick_markers(make_ref(setNames(as.numeric(!controls %in% drop), controls)), background), fold = f)
}))
w[, score_marker := gene %in% markers_any$gene]
fwrite(w[, .(gene, stage, season, class, log2FC_unadjusted, log2FC_adjusted = log2FC_adjusted_residual,
             fdr_unadjusted, fdr_adjusted = fdr_adjusted_residual, revealed_by_adjustment, sign_flip, score_marker)],
       file.path(OUT, "data/gene_delay_classification.csv"))

summ <- w[, .(n_unadjusted = sum(fdr_unadjusted < FDR), n_beyond_delay = sum(class == "beyond_delay"),
              n_explained_by_delay = sum(class == "explained_by_delay"), n_not_changed = sum(class == "not_changed"),
              n_revealed = sum(revealed_by_adjustment), n_sign_flip = sum(sign_flip)), by = .(season, stage)]
summ <- summ[order(season, match(stage, STAGES))]
summ[, pct_explained := 100 * n_explained_by_delay / n_unadjusted]
fwrite(summ, file.path(OUT, "data/adjusted_vs_unadjusted_summary.csv"))
print(summ); print(agree)

groups <- c(sHSP = "BRANCH:50", Photosystem_light_harvesting = "BRANCH:58",
            Cytosolic_translation = "BRANCH:55", Chloroplast_translation = "BRANCH:56")
g <- rbindlist(lapply(names(groups), function(n) {
  genes <- intersect(mem[set_id == groups[[n]]]$gene, rownames(E))
  w[gene %in% genes, .(group = n, n_genes = .N, beyond_delay = sum(class == "beyond_delay"),
                       explained_by_delay = sum(class == "explained_by_delay"), not_changed = sum(class == "not_changed"),
                       mean_log2FC_unadjusted = mean(log2FC_unadjusted), mean_log2FC_adjusted = mean(log2FC_adjusted_residual)),
    by = .(season, stage)]
}))
fwrite(g[order(group, season, match(stage, STAGES))], file.path(OUT, "data/headline_group_classes.csv"))
print(g[stage == "Harvest"])
cat("T4 COMPLETE\n")
