#!/usr/bin/env Rscript
# T7 part 1. Random-effects check and variance components.
ROOT_INIT <- Sys.getenv("GRBV_ROOT", if (dir.exists("Current_Analysis")) normalizePath(".") else normalizePath("../.."))
.libPaths(c(file.path(ROOT_INIT, "Reference_Resources/Rlibs_phase2"), .libPaths()))
source(file.path(ROOT_INIT, "Current_Analysis/simplified_reanalysis_2026-11/00_common.R"))
suppressPackageStartupMessages({library(BiocParallel); library(limma); library(variancePartition); library(ggplot2)})
source(file.path(ROOT, "Current_Analysis/simplified_reanalysis_2026-11/07_shared.R"))

vobj <- readRDS(file.path(OUT, "data/vobj_dream.rds"))
audit <- fread(file.path(OUT, "data/sample_audit.csv"))[match(meta$Sample, as.character(Sample))]
m$row <- factor(audit$row_number)
m$vine_season <- factor(paste(m$vine_id, m$season))
cat("Vines per row:\n"); print(unique(as.data.table(m)[, .(vine_id, row, infection)])[, .N, by = .(row, infection)][order(row, infection)])

models <- list(primary = ~ 0 + condition + (1 | vine_id),
               vine_season = ~ 0 + condition + (1 | vine_id) + (1 | vine_season),
               vine_season_row = ~ 0 + condition + row + (1 | vine_id) + (1 | vine_season))
ct <- fread(file.path(PRIOR, "Era6_contrast_definitions.csv"))
contrast_list <- setNames(ct$expression, ct$contrast)
res <- list(); cells <- list()
for (n in names(models)) {
  t0 <- Sys.time()
  L <- makeContrastsDream(models[[n]], m, contrasts = contrast_list)
  fit <- variancePartition::eBayes(dream(vobj, models[[n]], m, L = L, ddf = "Satterthwaite", BPPARAM = BP, quiet = TRUE))
  res[[n]] <- joint_counts(fit)
  cat(n, "seconds:", round(difftime(Sys.time(), t0, units = "secs")), "failures:", length(fit$errors), "\n"); print(res[[n]]$summary)
  cells[[n]] <- rbindlist(lapply(grep("^infection__", ct$contrast, value = TRUE), function(cn) {
    tt <- variancePartition::topTable(fit, coef = cn, number = Inf, sort.by = "none")
    data.table(gene = rownames(tt), contrast = cn, t = tt$t, fdr = tt$adj.P.Val)
  }))
}
out <- rbindlist(lapply(names(res), function(n) res[[n]]$summary[, model := n]))
for (n in names(models)) {
  mg <- merge(cells$primary, cells[[n]], by = c("gene", "contrast"), suffixes = c("_p", "_m"))
  a <- mg[, .(test = contrast[1], spearman_t = cor(t_p, t_m, method = "spearman"), n_fdr05 = sum(fdr_m < .05)), by = contrast]
  out <- rbind(out, a[, .(test, n_fdr05, spearman_t, model = n)], fill = TRUE)
}
fwrite(out, file.path(OUT, "data/dream_reffect_comparison.csv"))

form_vp <- ~ (1 | stage) + (1 | season) + (1 | infection) + (1 | vine_id) + (1 | vine_season) + (1 | row)
t0 <- Sys.time()
vp <- fitExtractVarPartModel(vobj, form_vp, m, BPPARAM = BP, quiet = TRUE)
cat("Variance partition seconds:", round(difftime(Sys.time(), t0, units = "secs")), "\n")
vd <- as.data.table(as.data.frame(vp), keep.rownames = "gene")
fwrite(vd, file.path(OUT, "data/variance_components.csv"))
long <- melt(vd, id.vars = "gene", variable.name = "component", value.name = "fraction")
print(long[, .(mean = mean(fraction), median = median(fraction)), by = component][order(-median)])
p <- ggplot(long, aes(reorder(component, -fraction, median), fraction)) +
  geom_violin(fill = "#9FB8B5", colour = NA, scale = "width") + geom_boxplot(width = .08, outlier.shape = NA, colour = "#7A1F3D") +
  labs(x = NULL, y = "Fraction of variance per gene") + theme_minimal(base_size = 9)
ggsave(file.path(OUT, "figures/SI_variance_components.png"), p, width = 7, height = 3.8, dpi = 300)
cat("T7 part 1 COMPLETE\n")
