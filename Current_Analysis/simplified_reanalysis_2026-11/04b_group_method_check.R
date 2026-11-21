#!/usr/bin/env Rscript
# T4 follow-up. The two adjustment methods disagree most at Harvest 2021, so the
# headline groups are checked under both, and gene classes are marked robust only
# when both methods call the gene.
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("GRBV_ROOT", if (dir.exists("Current_Analysis")) normalizePath(".") else normalizePath("../.."))
OUT <- file.path(ROOT, "Current_Analysis/simplified_reanalysis_2026-11")
mem <- fread(file.path(ROOT, "Current_Analysis/prior_era_inputs/Era7_pathway_membership.csv"))
a <- fread(file.path(OUT, "data/adjusted_dream_cell_contrasts.csv"))
w <- dcast(a, gene + stage + season ~ model, value.var = c("log2FC", "fdr"))
groups <- c(sHSP = "BRANCH:50", Photosystem_light_harvesting = "BRANCH:58",
            Cytosolic_translation = "BRANCH:55", Chloroplast_translation = "BRANCH:56")
r <- rbindlist(lapply(names(groups), function(n) {
  g <- unique(mem[set_id == groups[[n]]]$gene)
  w[gene %in% g, .(group = n, n_genes = .N,
    mean_residual = mean(log2FC_adjusted_residual), mean_covariate = mean(log2FC_adjusted_covariate),
    beyond_residual = sum(fdr_adjusted_residual < .05), beyond_covariate = sum(fdr_adjusted_covariate < .05),
    beyond_both = sum(fdr_adjusted_residual < .05 & fdr_adjusted_covariate < .05)), by = .(season, stage)]
}))
fwrite(r, file.path(OUT, "data/headline_group_method_agreement.csv"))
print(r[stage == "Harvest"])
s <- w[, .(beyond_residual = sum(fdr_adjusted_residual < .05), beyond_both = sum(fdr_adjusted_residual < .05 & fdr_adjusted_covariate < .05)), by = .(season, stage)]
fwrite(s, file.path(OUT, "data/beyond_delay_robust_counts.csv")); print(s)
