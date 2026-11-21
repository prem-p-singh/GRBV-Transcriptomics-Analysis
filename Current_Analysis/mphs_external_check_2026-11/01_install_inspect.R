# Step 1: install (done by R CMD INSTALL from vendor/MPhS at a pinned commit) and inspect MPhS.
# Run from Current_Analysis/mphs_external_check_2026-11.
.libPaths(c("Rlib", .libPaths()))
suppressPackageStartupMessages({library(MPhS); library(dplyr); library(tidyr)})
set.seed(20261107)
commit <- system("git -C vendor/MPhS rev-parse HEAD", intern = TRUE)
load(system.file("data", "MPhSdata.rda", package = "MPhS"))

# Smoke test: the README example (Sangiovese, Passerina; 4 stages x 3 replicates, replicate means)
data("RPKMdata", package = "MPhS")
exp_cond <- names(RPKMdata)[-1]; genes <- RPKMdata$gene_id
dts_vars <- data.frame(exp_cond) %>% separate(exp_cond, into = c("Cultivar", "Stage", "Replicate"), sep = "_")
dts <- as.data.frame(t(RPKMdata[, -1])); names(dts) <- genes; dts <- cbind(dts, dts_vars)
dts_means <- dts %>% group_by(Cultivar, Stage) %>% summarize(across(all_of(genes), mean), .groups = "drop")
out <- MPhStimepoints(data = dts_means, strata_var = "Cultivar", stage_var = "Stage")
smoke <- out$PCscores[, c("Cultivar", "Stage", "timepoint")]
smoke$Stage <- factor(smoke$Stage, c("Pea", "Touch", "Soft", "Harv"))
smoke <- smoke[order(smoke$Cultivar, smoke$Stage), ]
write.csv(smoke, "data/mphs_smoke_test.csv", row.names = FALSE)
mono <- tapply(smoke$timepoint, smoke$Cultivar, function(x) all(diff(x) > 0))

fmt <- function(x) paste(capture.output(x), collapse = "\n")
writeLines(c(
  "# Step 1: MPhS package inspection", "",
  sprintf("- Package: MPhS %s, GitHub sndmrc/MPhS, commit `%s`, installed into `Rlib` with `R CMD INSTALL` from `vendor/MPhS`.", packageVersion("MPhS"), commit),
  sprintf("- R: %s. Imports: dplyr %s, ggplot2 %s, ggrepel %s, rlang %s.", getRversion(), packageVersion("dplyr"), packageVersion("ggplot2"), packageVersion("ggrepel"), packageVersion("rlang")),
  "- Exported functions: `MPhStimepoints(data, strata_var, stage_var, scaling_type = \"scale\", geneID = \"VIT\")`, `MPhSscores(data, scaling_type, geneID)`, plus S3 `plot.MPhStimepoints`.",
  "",
  "## What the projection does (read from R/MPhStimepoints.R and R/MPhSscores.R)",
  "1. Input is a data frame with samples as rows and genes as columns named by V1 `VIT_` IDs (or V3 `Vitvi` IDs with `geneID = \"Vitvi\"`, converted through the bundled `v1v3` table). Extra columns hold strata and stage. Works per row: no replicate grouping is required (the README averages replicates, the function does not).",
  "2. Genes with any NA or with constant values are dropped; only genes in the 10,129-gene core set are kept.",
  "3. Scaling (`scaling_type`): `\"scale\"` (default) z-scores each gene **across the input samples**; `\"means_SDs\"` uses the reference means/SDs of the original dataset; `\"none\"` does nothing. No library-size normalization or log transform is done inside the package.",
  "4. Scaled data are multiplied by the stored PCA rotation (rows subset to the genes present). PCs 1, 2 and 5 are used; each is mapped through a linear fit from the raw to the smoothed reference scores.",
  "5. Each sample gets the **nearest of 30 fixed scale points** (`MPhSpts`, timepoint 1 to 30), by squared Euclidean distance. The output is therefore an integer from 1 to 30.",
  "",
  "## Bundled reference",
  sprintf("- Core genes: %d, ID format V1 `VIT_`. Examples: %s.", length(MPhSdata$means), paste(head(names(MPhSdata$means), 5), collapse = ", ")),
  sprintf("- Reference gene means: median %.2f, max %.0f; SDs median %.2f. These look like linear (not log) RPKM/FPKM values, matching the RPKM example data.", median(MPhSdata$means), max(MPhSdata$means), median(MPhSdata$sds)),
  sprintf("- PCA object class: %s; rotation %d x %d. Raw and smoothed reference scores: %d samples x %d PCs.", paste(class(MPhSdata$pca), collapse = "/"), nrow(MPhSdata$pca$rotation), ncol(MPhSdata$pca$rotation), nrow(MPhSdata$raw_data_scores), ncol(MPhSdata$raw_data_scores)),
  sprintf("- Scale points: %d (timepoints %d to %d). `v1v3` table: %d rows.", nrow(MPhSdata$MPhSpts), min(MPhSdata$MPhSpts$timepoint), max(MPhSdata$MPhSpts$timepoint), nrow(MPhSdata$v1v3)),
  "",
  "## Smoke test (README example, replicate means)",
  "```", fmt(smoke), "```",
  sprintf("Stages increase monotonically within each cultivar: %s. The README shows only a plot, not numbers, so this is checked for ordering, not for exact values.", paste(names(mono), mono, sep = " = ", collapse = "; ")),
  "",
  "## Answers to the plan's open questions",
  "- Gene IDs: V1 `VIT_` (12X); `Vitvi` V3 also accepted.",
  "- Input scale: linear expression values (reference is RPKM-like). The package does no normalization of its own.",
  "- Internal normalization: only per-gene scaling (default z-score across the input samples).",
  "- Output: per-sample integer stage 1 to 30; no replicate groups needed.",
  "",
  "## Concerns logged for the author (not acted on)",
  "- **Default scaling is relative to the input set.** With `\"scale\"`, each gene is centred on the mean of whatever libraries are passed in, so projecting the 63 healthy libraries alone and the 327 libraries together can give different stages for the same library. Proposal: in all steps, project all 327 libraries in one call so every library shares one scaling, and run the healthy gate on the healthy rows of that call. Projecting healthy libraries alone is reported as a sensitivity check.",
  "- **The output is a whole number from 1 to 30.** A delay of 0.15 in-house stage units may be smaller than one MPhS step, so a small 2022 delay could round to zero for most libraries. This limits the precision, not the direction, of the check.",
  "- **3' TAG-seq vs full-length RPKM.** The reference used full-length RNA-seq. Gene-wise z-scoring (default) removes per-gene level differences such as length effects, which is why the plan's no-length-correction CPM is a reasonable input."
), "logs/01_package_inspection.md")
writeLines(capture.output(sessionInfo()), "logs/sessionInfo.txt")
cat("monotonic:", mono, "\n"); print(smoke)
