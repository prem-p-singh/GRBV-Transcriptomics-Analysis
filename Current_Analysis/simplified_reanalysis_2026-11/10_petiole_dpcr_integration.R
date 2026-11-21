#!/usr/bin/env Rscript
# 10_petiole_dpcr_integration.R
# Cross-tissue integration of companion leaf-petiole QuantStudio 3D dPCR viral DNA
# (Singh et al., 2026, Arch Virol 171:183; Raw_Data/Petiole_dPCR_Viral_Titer_VineID.xlsx)
# with berry 3' TAG-seq across the exact same 36 GRBV-positive Merlot vines (R{row}V{vine}).

ROOT <- Sys.getenv("GRBV_ROOT", if (dir.exists("Current_Analysis")) normalizePath(".") else normalizePath("../.."))
.libPaths(c(file.path(ROOT, "Reference_Resources/Rlibs_phase2"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(variancePartition)
  library(BiocParallel)
})

SIMP <- file.path(ROOT, "Current_Analysis/simplified_reanalysis_2026-11/data")
TAB  <- file.path(ROOT, "Manuscript/tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

STAGES <- c("Pre_veraison", "Veraison", "Post_veraison", "Harvest")

dpcr <- fread(file.path(SIMP, "petiole_dpcr_vine_stage_season.csv"))[!is.na(adjusted_copy_no)]
dpcr[, `:=`(
  season = as.character(season),
  stage = factor(stage, levels = STAGES),
  log10_petiole = log10(adjusted_copy_no),
  log2_petiole = log2(adjusted_copy_no)
)]

# 1. Stage x season summary table
p_sum <- dpcr[, .(
  n_chips = .N,
  mean_log10 = mean(log10_petiole),
  sd_log10 = sd(log10_petiole),
  se_log10 = sd(log10_petiole) / sqrt(.N),
  geomean_copies = 10^mean(log10_petiole),
  median_copies = as.numeric(median(adjusted_copy_no)),
  iqr_low = as.numeric(quantile(adjusted_copy_no, 0.25)),
  iqr_high = as.numeric(quantile(adjusted_copy_no, 0.75))
), by = .(season, stage)][order(season, stage)]
fwrite(p_sum, file.path(SIMP, "petiole_dpcr_stage_season_summary.csv"))
fwrite(p_sum, file.path(TAB, "petiole_dpcr_stage_season_summary.csv"))

# 2. Merge with 264 infected berry RNA-seq libraries by exact R{row}V{vine}
sa   <- fread(file.path(SIMP, "sample_audit.csv"))[, season := as.character(season)]
wlib <- fread(file.path(SIMP, "within_infected_library_values.csv"))[, season := as.character(season)]
vobj <- readRDS(file.path(SIMP, "vobj_dream.rds"))

viral_ids <- c("N761_gp1", "N761_gp2", "N761_gp3", "N761_gp4", "N761_gp5")
v_expr <- as.data.table(t(vobj$E[viral_ids, , drop = FALSE]), keep.rownames = "Sample")
setnames(v_expr, viral_ids, c("V2", "V1", "V3", "C2", "C1"))
v_expr[, `:=`(C2_minus_V2 = C2 - V2, C2_minus_C1 = C2 - C1)]

m <- merge(wlib, sa[, .(Sample, viral_cpm, library_size)], by = "Sample")
m <- merge(m, v_expr, by = "Sample")
m <- merge(
  m,
  dpcr[, .(vine = vine_id, season, stage = as.character(stage),
           adjusted_copy_no, copies_per_ul, log10_petiole, log2_petiole)],
  by = c("vine", "season", "stage"),
  all.x = TRUE
)
m[, stage := factor(stage, levels = STAGES)]
fwrite(m, file.path(SIMP, "paired_petiole_berry_libraries.csv"))

# Update titer_vs_viralcpm.csv with verified Vine-ID Spearman & Pearson correlations
mp <- m[!is.na(log10_petiole)]
cor_rows <- rbind(
  data.table(
    level = "overall", season = "all", stage = "all", n = nrow(mp),
    spearman_rho = cor(mp$adjusted_copy_no, mp$viral_cpm, method = "spearman"),
    pearson_r_log10 = cor(mp$log10_petiole, mp$log10_viral),
    p_pearson = cor.test(mp$log10_petiole, mp$log10_viral)$p.value
  ),
  mp[, .(
    level = "cell", n = .N,
    spearman_rho = cor(adjusted_copy_no, viral_cpm, method = "spearman"),
    pearson_r_log10 = cor(log10_petiole, log10_viral),
    p_pearson = cor.test(log10_petiole, log10_viral)$p.value
  ), by = .(season, stage)][order(season, stage)]
)
fwrite(cor_rows, file.path(SIMP, "titer_vs_viralcpm.csv"))

# 3. Cross-tissue mixed-model summary tests
m_anova <- anova(lmer(log10_petiole ~ stage * season + (1 | vine_id), data = dpcr))
m21 <- summary(lmer(log10_petiole ~ stage + (1 | vine_id), data = dpcr[season == "2021"]))$coefficients
m22 <- summary(lmer(log10_petiole ~ stage + (1 | vine_id), data = dpcr[season == "2022"]))$coefficients
m_pre <- summary(lmer(log10_petiole ~ season + (1 | vine_id), data = dpcr[stage == "Pre_veraison"]))$coefficients
m_har <- summary(lmer(log10_petiole ~ season + (1 | vine_id), data = dpcr[stage == "Harvest"]))$coefficients

dpcr_w <- dcast(dpcr, vine_id + season ~ stage, value.var = "log10_petiole")
dpcr_w[, `:=`(
  pet_early = rowMeans(cbind(Pre_veraison, Veraison), na.rm = TRUE),
  pet_mean  = rowMeans(cbind(Pre_veraison, Veraison, Post_veraison, Harvest), na.rm = TRUE)
)]
har <- merge(m[stage == "Harvest"], dpcr_w, by.x = c("vine", "season"), by.y = c("vine_id", "season"))

fit_har_mean  <- summary(lmer(developmental_offset ~ season + pet_mean + (1 | vine), data = har))$coefficients["pet_mean", ]
fit_har_early <- summary(lmer(developmental_offset ~ season + pet_early + (1 | vine), data = har))$coefficients["pet_early", ]
har21 <- har[season == "2021" & !is.na(Harvest)]
har22 <- har[season == "2022" & !is.na(Harvest)]
r21 <- cor.test(har21$Harvest, har21$developmental_offset)
rho21 <- cor.test(har21$Harvest, har21$developmental_offset, method = "spearman", exact = FALSE)
rho22 <- cor.test(har22$Harvest, har22$developmental_offset, method = "spearman", exact = FALSE)

fit_berry_couple <- summary(lmer(log10_viral ~ stage * season + log10_petiole + (1 | vine), data = mp))$coefficients["log10_petiole", ]
fit_c2_couple    <- summary(lmer(C2 ~ stage * season + log10_petiole + (1 | vine), data = mp))$coefficients["log10_petiole", ]
fit_shsp_couple  <- summary(lmer(shsp_adjusted_residual ~ stage * season + log10_petiole + (1 | vine), data = mp))$coefficients["log10_petiole", ]

# 4. Joint genome-wide dream scan (petiole DNA vs berry viral RNA)
host_genes <- grep("^N761_", rownames(vobj$E), value = TRUE, invert = TRUE)
vobj_sub <- vobj[host_genes, match(mp$Sample, colnames(vobj$E))]
meta_sub <- as.data.frame(mp); rownames(meta_sub) <- meta_sub$Sample
param <- MulticoreParam(workers = 4)

fit_pet_only <- eBayes(dream(vobj_sub, ~ stage * season + log10_petiole + (1 | vine), meta_sub, BPPARAM = param))
tt_pet_only  <- as.data.table(topTable(fit_pet_only, coef = "log10_petiole", number = Inf), keep.rownames = "gene")

fit_joint <- eBayes(dream(vobj_sub, ~ stage * season + log10_petiole + log10_viral + (1 | vine), meta_sub, BPPARAM = param))
tt_j_pet  <- as.data.table(topTable(fit_joint, coef = "log10_petiole", number = Inf), keep.rownames = "gene")
tt_j_ber  <- as.data.table(topTable(fit_joint, coef = "log10_viral", number = Inf), keep.rownames = "gene")

joint_dt <- merge(
  tt_j_pet[, .(gene, beta_petiole = logFC, t_petiole = t, p_petiole = P.Value, fdr_petiole = adj.P.Val)],
  tt_j_ber[, .(gene, beta_berry_viral = logFC, t_berry_viral = t, p_berry_viral = P.Value, fdr_berry_viral = adj.P.Val)],
  by = "gene"
)
fwrite(joint_dt, file.path(SIMP, "joint_petiole_berry_host_dream.csv"))

tests_dt <- data.table(
  test_id = c(
    "petiole_stage_x_season_F",
    "petiole_2021_harvest_vs_pre_log10",
    "petiole_2022_harvest_vs_pre_log10",
    "petiole_pre_veraison_2022_vs_2021_log10",
    "petiole_harvest_2022_vs_2021_log10",
    "berry_viral_vs_petiole_RE",
    "berry_C2_vs_petiole_RE",
    "harvest_dev_offset_vs_season_mean_petiole",
    "harvest_dev_offset_vs_early_petiole",
    "harvest_2021_dev_offset_vs_harvest_petiole_pearson",
    "harvest_2021_dev_offset_vs_harvest_petiole_spearman",
    "harvest_2022_dev_offset_vs_harvest_petiole_spearman",
    "berry_shsp_residual_vs_petiole_RE",
    "host_genes_petiole_only_fdr05",
    "host_genes_joint_petiole_fdr05",
    "host_genes_joint_berry_viral_fdr05"
  ),
  estimate = c(
    m_anova["stage:season", "F value"],
    m21["stageHarvest", "Estimate"],
    m22["stageHarvest", "Estimate"],
    m_pre["season2022", "Estimate"],
    m_har["season2022", "Estimate"],
    fit_berry_couple["Estimate"],
    fit_c2_couple["Estimate"],
    fit_har_mean["Estimate"],
    fit_har_early["Estimate"],
    r21$estimate,
    rho21$estimate,
    rho22$estimate,
    fit_shsp_couple["Estimate"],
    sum(tt_pet_only$adj.P.Val < 0.05),
    sum(tt_j_pet$adj.P.Val < 0.05),
    sum(tt_j_ber$adj.P.Val < 0.05)
  ),
  p_or_min_fdr = c(
    m_anova["stage:season", "Pr(>F)"],
    m21["stageHarvest", "Pr(>|t|)"],
    m22["stageHarvest", "Pr(>|t|)"],
    m_pre["season2022", "Pr(>|t|)"],
    m_har["season2022", "Pr(>|t|)"],
    fit_berry_couple["Pr(>|t|)"],
    fit_c2_couple["Pr(>|t|)"],
    fit_har_mean["Pr(>|t|)"],
    fit_har_early["Pr(>|t|)"],
    r21$p.value,
    rho21$p.value,
    rho22$p.value,
    fit_shsp_couple["Pr(>|t|)"],
    min(tt_pet_only$adj.P.Val),
    min(tt_j_pet$adj.P.Val),
    min(tt_j_ber$adj.P.Val)
  )
)
fwrite(tests_dt, file.path(SIMP, "petiole_berry_cross_tissue_tests.csv"))
fwrite(tests_dt, file.path(TAB, "petiole_berry_cross_tissue_tests.csv"))
cat("Saved all cross-tissue petiole + berry integration tables successfully!\n")
print(tests_dt)
