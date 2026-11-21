#!/usr/bin/env Rscript
# Step 5: does the sHSP deficit beyond the delay survive when the delay is measured on MPhS?
# Same voom object, design and contrasts as simplified_reanalysis_2026-11/04_adjusted_dream.R
# (its "adjusted_covariate" model), with the MPhS offset in place of the in-house offset.
# Run from the repository root: Rscript Current_Analysis/mphs_external_check_2026-11/05_mphs_adjusted_de.R
ROOT_INIT <- Sys.getenv("GRBV_ROOT", if (dir.exists("Current_Analysis")) normalizePath(".") else normalizePath("../.."))
.libPaths(c(file.path(ROOT_INIT, "Reference_Resources/Rlibs_phase2"), .libPaths()))
source(file.path(ROOT_INIT, "Current_Analysis/simplified_reanalysis_2026-11/00_common.R"))   # read-only: meta, E, mem, STAGES, SEASONS, CORES
suppressPackageStartupMessages({library(BiocParallel); library(edgeR); library(limma); library(variancePartition)})
HERE <- file.path(ROOT_INIT, "Current_Analysis/mphs_external_check_2026-11")
STAGE_CODES <- c(Pre_veraison = "Pre", Veraison = "Veraison", Post_veraison = "Post", Harvest = "Harvest")
FDR <- 0.05
BP <- MulticoreParam(CORES, progressbar = FALSE)

m <- as.data.frame(meta); rownames(m) <- m$Sample
m$condition <- factor(paste(m$infection, STAGE_CODES[as.character(m$Maturity_stage)], m$season, sep = "_"))
m$vine_id <- factor(m$vine_id)
m$stage <- factor(as.character(m$Maturity_stage), levels = STAGES)

# MPhS offset: library stage minus the mean healthy MPhS stage of its nominal stage x season.
mp <- fread(file.path(HERE, "data/mphs_all_projection.csv"), colClasses = c(Sample = "character", season = "character"))
mp <- mp[match(m$Sample, Sample)]
stopifnot(!anyNA(mp$mphs), identical(mp$Sample, m$Sample), identical(mp$stage, as.character(m$stage)))
hmean <- mp[infection == "RB(-)", .(h = mean(mphs)), by = .(season, stage)]
m$mphs_offset <- mp$mphs - hmean$h[match(paste(mp$season, mp$stage), paste(hmean$season, hmean$stage))]
stopifnot(!anyNA(m$mphs_offset))

vobj <- readRDS(file.path(OUT, "data/vobj_dream.rds"))
stopifnot(identical(colnames(vobj$E), m$Sample), max(abs(vobj$E - E)) < 1e-8)

contrasts_all <- c()
for (y in SEASONS) for (s in STAGES) contrasts_all[paste("infection", STAGE_CODES[[s]], y, sep = "__")] <-
  sprintf("conditionpos_%s_%s - conditionneg_%s_%s", STAGE_CODES[[s]], y, STAGE_CODES[[s]], y)
form <- ~ 0 + condition + mphs_offset:stage + (1 | vine_id)
L <- makeContrastsDream(form, m, contrasts = contrasts_all)

t0 <- Sys.time()
fit <- dream(vobj, form, m, L = L, ddf = "Satterthwaite", BPPARAM = BP, quiet = TRUE)
fit <- variancePartition::eBayes(fit)
cat("dream seconds:", round(difftime(Sys.time(), t0, units = "secs")), "failures:", length(fit$errors), "\n")
tab <- rbindlist(lapply(names(contrasts_all), function(cn) {
  tt <- variancePartition::topTable(fit, coef = cn, number = Inf, sort.by = "none", adjust.method = "BH")
  f <- strsplit(cn, "__", fixed = TRUE)[[1]]
  data.table(gene = rownames(tt), stage = names(STAGE_CODES)[match(f[2], STAGE_CODES)], season = f[3],
             log2FC = tt$logFC, t = tt$t, fdr = tt$adj.P.Val)
}))
stopifnot(nrow(tab) == 14604 * 8)
fwrite(tab, file.path(HERE, "data/mphs_adjusted_cell_contrasts.csv.gz"))

# In-house results for comparison: covariate model (like for like) and residual model (headline file).
inh <- fread(file.path(OUT, "data/adjusted_dream_cell_contrasts.csv"))[model == "adjusted_covariate"]
inh[, season := as.character(season)]
unadj <- fread(file.path(OUT, "data/unadjusted_dream_cell_contrasts.csv"))[, season := as.character(season)]
head_inh <- fread(file.path(OUT, "data/headline_group_classes.csv"))[, season := as.character(season)]

groups <- c(sHSP = "BRANCH:50", Photosystem_light_harvesting = "BRANCH:58",
            Cytosolic_translation = "BRANCH:55", Chloroplast_translation = "BRANCH:56")
summ <- rbindlist(lapply(names(groups), function(n) {
  genes <- intersect(mem[set_id == groups[[n]]]$gene, rownames(E))
  a <- tab[gene %in% genes & stage == "Harvest", .(n_genes = .N, mphs_mean_log2FC = mean(log2FC), mphs_n_fdr05 = sum(fdr < FDR),
                                                   mphs_n_down = sum(fdr < FDR & log2FC < 0), mphs_n_up = sum(fdr < FDR & log2FC > 0)), by = season]
  b <- inh[gene %in% genes & stage == "Harvest", .(inhouse_cov_mean_log2FC = mean(log2FC), inhouse_cov_n_fdr05 = sum(fdr < FDR)), by = season]
  u <- unadj[gene %in% genes & stage == "Harvest", .(unadjusted_mean_log2FC = mean(log2FC), unadjusted_n_fdr05 = sum(fdr < FDR)), by = season]
  h <- head_inh[group == n & stage == "Harvest", .(season, inhouse_resid_mean_log2FC = mean_log2FC_adjusted, inhouse_resid_n_beyond_delay = beyond_delay)]
  r <- Reduce(function(x, y) merge(x, y, by = "season"), list(a, b, u, h)); r[, group := n]; r
}))
setcolorder(summ, c("group", "season"))
fwrite(summ, file.path(HERE, "data/mphs_adjusted_headline_groups.csv"))

# Genome-wide agreement with the in-house covariate model, per cell
j <- merge(tab, inh[, .(gene, stage, season, i_lfc = log2FC, i_t = t, i_fdr = fdr)], by = c("gene", "stage", "season"))
ag <- j[, .(spearman_t = cor(t, i_t, method = "spearman"), pearson_log2FC = cor(log2FC, i_lfc),
            n_mphs = sum(fdr < FDR), n_inhouse = sum(i_fdr < FDR), n_both = sum(fdr < FDR & i_fdr < FDR)), by = .(season, stage)]
ag[, jaccard := n_both / (n_mphs + n_inhouse - n_both)]
fwrite(ag[order(season, match(stage, STAGES))], file.path(HERE, "data/mphs_vs_inhouse_adjusted_agreement.csv"))
writeLines(capture.output(sessionInfo()), file.path(HERE, "logs/sessionInfo_step5.txt"))
print(summ, digits = 3); print(ag, digits = 3)
cat("STEP 5 COMPLETE\n")
