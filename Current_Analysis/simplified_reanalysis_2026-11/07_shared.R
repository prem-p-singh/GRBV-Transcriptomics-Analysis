# Model setup and joint-test counts shared by the T7 scripts.
STAGE_CODES <- c(Pre_veraison = "Pre", Veraison = "Veraison", Post_veraison = "Post", Harvest = "Harvest")
BP <- MulticoreParam(CORES, progressbar = FALSE)
m <- as.data.frame(meta)
rownames(m) <- m$Sample
make_condition <- function(g) factor(paste(g, STAGE_CODES[as.character(m$Maturity_stage)], m$season, sep = "_"))
m$condition <- make_condition(m$infection)
m$vine_id <- factor(m$vine_id)
m$stage <- factor(as.character(m$Maturity_stage), levels = STAGES)
m$vine_season <- factor(paste(m$vine_id, m$season))

JOINT <- list(
  infection_stage = paste0("infection_stage_average__", c("Veraison", "Post", "Harvest")),
  infection_stage_season = paste0("three_way__", c("Veraison", "Post", "Harvest")))

# Same F test as 22_longitudinal_dream.R. only_season_specific counts genes found
# by the three-way test and not by the season-averaged test.
joint_counts <- function(fit) {
  fdr <- lapply(JOINT, function(cf) {
    tt <- variancePartition::topTable(fit, coef = cf, number = Inf, sort.by = "none", adjust.method = "none")
    p.adjust(tt$P.Value, "BH")
  })
  list(summary = data.table(test = c("infection_stage", "infection_stage_season", "only_season_specific"),
                            n_fdr05 = c(sum(fdr[[1]] < .05, na.rm = TRUE), sum(fdr[[2]] < .05, na.rm = TRUE),
                                        sum(fdr[[2]] < .05 & fdr[[1]] >= .05, na.rm = TRUE))), fdr = fdr)
}
