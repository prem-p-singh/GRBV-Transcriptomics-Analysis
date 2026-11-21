#!/usr/bin/env Rscript
# Bundles what the FARM permutation needs: voom object, model frame, contrasts.
ROOT_INIT <- Sys.getenv("GRBV_ROOT", if (dir.exists("Current_Analysis")) normalizePath(".") else normalizePath("../../.."))
.libPaths(c(file.path(ROOT_INIT, "Reference_Resources/Rlibs_phase2"), .libPaths()))
source(file.path(ROOT_INIT, "Current_Analysis/simplified_reanalysis_2026-11/00_common.R"))
suppressPackageStartupMessages({library(BiocParallel); library(limma); library(variancePartition)})
source(file.path(ROOT, "Current_Analysis/simplified_reanalysis_2026-11/07_shared.R"))
ct <- fread(file.path(PRIOR, "Era6_contrast_definitions.csv"))
vobj <- readRDS(file.path(OUT, "data/vobj_dream.rds"))
keep <- c("Sample", "season", "Maturity_stage", "infection", "vine_id", "condition", "stage", "vine_season")
saveRDS(list(vobj = vobj, m = m[, keep], contrasts = setNames(ct$expression, ct$contrast), joint = JOINT, stage_codes = STAGE_CODES, seed = SEED),
        file.path(OUT, "farm/perm_inputs.rds"))
cat("Saved inputs\n")
