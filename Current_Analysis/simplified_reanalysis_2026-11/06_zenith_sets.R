#!/usr/bin/env Rscript
# T6. One competitive gene-set test from the dream fits: zenith with
# inter.gene.cor = NA, which estimates the correlation per set from residuals.
ROOTDIR <- Sys.getenv("GRBV_ROOT", if (dir.exists("Current_Analysis")) normalizePath(".") else normalizePath("../.."))
.libPaths(c(file.path(ROOTDIR, "Current_Analysis/simplified_reanalysis_2026-11/Rlib"), file.path(ROOTDIR, "Reference_Resources/Rlibs_phase2"), .libPaths()))
source(file.path(ROOTDIR, "Current_Analysis/simplified_reanalysis_2026-11/00_common.R"))
suppressPackageStartupMessages({library(parallel); library(variancePartition); library(zenith)})
fits <- readRDS(file.path(OUT, "data/dream_fits.rds"))
stopifnot(!is.null(fits$primary$residuals), !is.null(fits$adjusted$residuals))
cov <- fread(file.path(PRIOR, "Era7_pathway_coverage.csv"))
eligible <- cov[eligible == TRUE]$set_id
m2 <- mem[set_id %in% eligible & gene %in% rownames(E), .(gene = unique(gene)), by = set_id]
index <- split(m2$gene, m2$set_id)
index <- index[lengths(index) >= 5 & lengths(index) <= 500]
n_branch <- sum(startsWith(names(index), "BRANCH:"))
cat("Sets tested:", length(index), " branches:", n_branch, "\n")
stopifnot(length(index) == 1369, n_branch == 59)

run <- function(fit, label) {
  coefs <- grep("^infection__", colnames(coef(fit)), value = TRUE)
  stopifnot(length(coefs) == 8)
  r <- rbindlist(mclapply(coefs, function(cn) {
    z <- zenith(fit, cn, index, inter.gene.cor = NA, progressbar = FALSE)
    z <- as.data.table(z, keep.rownames = "set_id")
    f <- strsplit(cn, "__", fixed = TRUE)[[1]]
    z[, `:=`(stage_code = f[2], season = f[3])]
    z
  }, mc.cores = min(8L, CORES)))
  stopifnot(nrow(r) == 8 * length(index))
  r[, FDR_all := p.adjust(PValue, "BH")]
  r[startsWith(set_id, "BRANCH:"), FDR_branches := p.adjust(PValue, "BH")]  # prespecified family, 59 x 8
  r <- merge(r, cov[, .(set_id, name, collection)], by = "set_id")
  fwrite(r, file.path(OUT, sprintf("data/zenith_%s.csv", label)))
  r
}
un <- run(fits$primary, "unadjusted")
ad <- run(fits$adjusted, "adjusted")

# Agreement with the Era 7 fry results for the branches named in the manuscript.
named <- c("BRANCH:50", "BRANCH:55", "BRANCH:56", "BRANCH:58")
code <- c(Pre_veraison = "Pre", Veraison = "Veraison", Post_veraison = "Post", Harvest = "Harvest")
fry <- fread(file.path(PRIOR, "Era7_pathway_fry_all_cells.csv"))[startsWith(set_id, "BRANCH:")]
fry[, c("stage", "season") := tstrsplit(cell, "__", fixed = TRUE)]
fry[, stage_code := code[stage]]
keep <- function(z, tag) z[startsWith(set_id, "BRANCH:"),
  setNames(list(set_id, stage_code, season, Direction, PValue, FDR_all, FDR_branches, delta, Correlation),
           c("set_id", "stage_code", "season", paste0(tag, c("_direction", "_p", "_fdr", "_fdr_branches", "_delta", "_correlation"))))]
cmp <- merge(merge(keep(un, "zenith_unadj"), keep(ad, "zenith_adj"), by = c("set_id", "stage_code", "season")),
             fry[, .(set_id, stage_code, season, name, fry_direction = Direction, fry_fdr = FDR_all_cells, fry_fdr_mixed = FDR_Mixed_all_cells)],
             by = c("set_id", "stage_code", "season"))
cmp[, `:=`(same_direction = zenith_unadj_direction == fry_direction,
           both_sig = zenith_unadj_fdr < .05 & fry_fdr < .05)]
fwrite(cmp, file.path(OUT, "data/settest_old_vs_new.csv"))
cat("Branch-cell rows:", nrow(cmp), " same direction:", mean(cmp$same_direction), "\n")
cat("Branch family FDR<.05: unadjusted", cmp[, sum(zenith_unadj_fdr_branches < .05)], " adjusted", cmp[, sum(zenith_adj_fdr_branches < .05)], "\n")
cat("Sig FDR<.05: zenith unadjusted", cmp[, sum(zenith_unadj_fdr < .05)], " zenith adjusted", cmp[, sum(zenith_adj_fdr < .05)],
    " fry", cmp[, sum(fry_fdr < .05)], " both", sum(cmp$both_sig), "\n")
print(cmp[set_id %in% named & stage_code == "Harvest",
          .(name, season, zenith_unadj_p, zenith_unadj_fdr_branches, zenith_adj_p, zenith_adj_fdr_branches, fry_direction, fry_fdr)])
cat("T6 COMPLETE\n")
