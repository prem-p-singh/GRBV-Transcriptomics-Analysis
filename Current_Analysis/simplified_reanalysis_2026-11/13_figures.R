#!/usr/bin/env Rscript
# T13. Main figures. Every plotted value comes from a data/ CSV of this revision.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(scales); library(ragg)})
ROOT <- Sys.getenv("GRBV_ROOT", if (dir.exists("Current_Analysis")) normalizePath(".") else normalizePath("../.."))
OUT <- file.path(ROOT, "Current_Analysis/simplified_reanalysis_2026-11")
D <- file.path(OUT, "data"); FIG <- file.path(OUT, "figures")
C <- c("2021" = "#087F8C", "2022" = "#8C2846")   # project palette
SHAPE <- c("2021" = 16, "2022" = 17)
INK <- "#454B50"; GRID <- "#E9EDF0"; MID <- "#8A9299"
STAGES <- c("Pre_veraison", "Veraison", "Post_veraison", "Harvest")
SHORT <- c("Pre", "Ver", "Post", "Har")
BASE <- 8
theme_fig <- function() theme_classic(base_size = BASE, base_family = "Times New Roman") +
  theme(axis.line = element_line(colour = MID, linewidth = .3), axis.ticks = element_line(colour = MID, linewidth = .3),
        axis.text = element_text(colour = "black", size = 7), panel.grid.major.y = element_line(colour = GRID, linewidth = .3),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 7.5),
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 7), legend.title = element_text(size = 7),
        plot.title = element_text(size = 8.5, face = "bold"), plot.subtitle = element_text(size = 7, colour = INK),
        plot.tag = element_text(size = 11, face = "bold"), plot.margin = margin(3, 5, 3, 3))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 600, device = agg_png, bg = "white")
  ggsave(file.path(FIG, paste0(name, ".svg")), p, width = w, height = h, bg = "white")
  cat(name, "saved\n")
}
MINUS <- label_number(style_negative = "minus")
season_scales <- list(scale_colour_manual(values = C, name = "Season"), scale_shape_manual(values = SHAPE, name = "Season"))

## Figure 1
pred <- fread(file.path(D, "dev_score_predictions.csv"), colClasses = c(season = "character"))
hold <- fread(file.path(D, "dev_score_heldout_season.csv"), colClasses = c(season = "character"))
val <- fread(file.path(D, "dev_score_validation.csv"))
inf <- fread(file.path(D, "dev_displacement_inference.csv"), colClasses = c(season = "character"))
valid_panel <- function(d, label, title) {
  r <- val[validation == label & score == "new_score"]
  ggplot(d, aes(stage_number, new_score, colour = season, shape = season)) +
    geom_abline(linetype = "dashed", colour = MID, linewidth = .35) +
    geom_point(position = position_jitter(width = .12, height = 0, seed = 1), size = 1.1, alpha = .85) +
    annotate("text", x = -.2, y = 3.3, hjust = 0, vjust = 1, size = 2.5, family = "Times New Roman",
             label = sprintf("%d of %d nearest stage\nRMSE %.3f", r$nearest_correct, r$n, r$RMSE)) +
    season_scales + scale_x_continuous(breaks = 0:3, labels = SHORT) + scale_y_continuous(labels = MINUS) + coord_cartesian(ylim = c(-.25, 3.45)) +
    labs(title = title, x = "Recorded stage", y = "Developmental score") + theme_fig() + theme(legend.position = "none")
}
a <- valid_panel(pred[infection == "neg"], "held_out_vine", "Whole vine withheld")
b <- valid_panel(hold, "held_out_season", "Whole season withheld")
inf[, stage := factor(stage, STAGES)]
c1 <- ggplot(inf, aes(stage, estimate, colour = season, shape = season)) +
  geom_hline(yintercept = 0, colour = MID, linewidth = .35) +
  geom_errorbar(aes(ymin = nested_ci_low, ymax = nested_ci_high), width = 0, linewidth = .45, position = position_dodge(.45)) +
  geom_point(size = 1.9, position = position_dodge(.45)) +
  geom_text(aes(y = nested_ci_low - .05, label = ifelse(perm_fdr < .0005, "<0.0005", sprintf("%.3g", perm_fdr))), vjust = 1, size = 2.3,
            family = "Times New Roman", position = position_dodge(.45), show.legend = FALSE) +
  season_scales + scale_x_discrete(labels = c("Pre-veraison", "Veraison", "Post-veraison", "Harvest")) +
  scale_y_continuous(labels = MINUS) + coord_cartesian(ylim = c(-1.15, .25)) +
  labs(title = "Infected minus healthy score", subtitle = "Bars: nested bootstrap 95% CI. Numbers: permutation FDR.",
       x = NULL, y = "Difference (stage units)") + theme_fig()
save_fig((a | b) / (c1 + theme(legend.position = "bottom")) + plot_layout(heights = c(1, 1.05)) + plot_annotation(tag_levels = "A"), "Figure1", 7, 6)

## Figure 2
GROUPS <- c(sHSP = "BRANCH:50", Cytosolic_translation = "BRANCH:55", Chloroplast_translation = "BRANCH:56", Photosystem_light_harvesting = "BRANCH:58")
LAB <- c(sHSP = "Small heat-shock\nproteins", Cytosolic_translation = "Cytosolic\ntranslation",
         Chloroplast_translation = "Chloroplast\ntranslation", Photosystem_light_harvesting = "Photosystem and\nlight harvesting")
mem <- fread(file.path(ROOT, "Current_Analysis/prior_era_inputs/Era7_pathway_membership.csv"))
cl <- fread(file.path(D, "gene_delay_classification.csv"), colClasses = c(season = "character"))[stage == "Harvest"]
g <- rbindlist(lapply(names(GROUPS), function(n) cl[gene %in% mem[set_id == GROUPS[[n]]]$gene][, group := n]))
g[, group := factor(LAB[group], LAB)]
g[, class2 := factor(fifelse(class == "beyond_delay", "Beyond delay", "Explained or unchanged"), c("Beyond delay", "Explained or unchanged"))]
p2a <- ggplot(g, aes(log2FC_unadjusted, log2FC_adjusted, colour = class2)) +
  geom_abline(linetype = "dashed", colour = MID, linewidth = .3) + geom_hline(yintercept = 0, colour = GRID) + geom_vline(xintercept = 0, colour = GRID) +
  geom_point(size = .8, alpha = .8, shape = 16) +
  facet_grid(season ~ group, scales = "free") + scale_x_continuous(labels = MINUS) + scale_y_continuous(labels = MINUS) +
  scale_colour_manual(values = c("Beyond delay" = "#8C2846", "Explained or unchanged" = "#A4AEB0"), name = NULL) +
  labs(title = "Gene effects at harvest before and after removing the shared delay",
       subtitle = "Points on the dashed line are unchanged by adjustment; points nearer zero fit the delay.",
       x = "log2 fold change, unadjusted", y = "log2 fold change, adjusted") + theme_fig() + theme(legend.position = "bottom", panel.grid.major.y = element_blank())
nl <- fread(file.path(D, "matched_null_percentiles.csv"), colClasses = c(season = "character"))[stage == "Harvest"]
nl[, group := factor(LAB[group], rev(LAB))]
nl[, statistic := factor(fifelse(statistic == "mean_adjusted_effect", "Mean adjusted effect (log2)", "Score offset (stage units)"),
                         c("Mean adjusted effect (log2)", "Score offset (stage units)"))]
p2b <- ggplot(nl, aes(y = group, colour = season, shape = season)) +
  geom_vline(xintercept = 0, colour = MID, linewidth = .3) +
  geom_linerange(aes(xmin = null_q025, xmax = null_q975), linewidth = 2.2, alpha = .28, position = position_dodge(.6)) +
  geom_point(aes(x = real_value), size = 1.8, position = position_dodge(.6)) +
  facet_wrap(~ statistic, scales = "free_x") + season_scales + scale_x_continuous(labels = MINUS) +
  labs(title = "Real group against matched random sets", subtitle = "Band: central 95% of 1,999 matched sets. Point: real group. A point outside the band passes.",
       x = NULL, y = NULL) + theme_fig() + theme(legend.position = "bottom", panel.grid.major.y = element_blank())
save_fig(p2a / p2b + plot_layout(heights = c(1.5, 1)) + plot_annotation(tag_levels = "A"), "Figure2", 7, 7.2)

## Figure 3
ov <- fread(file.path(D, "interaction_overlap.csv"))[viral == FALSE]
bars <- data.table(set = factor(c("Infection × stage only", "Both tests", "Infection × stage × season only"),
                                rev(c("Infection × stage only", "Both tests", "Infection × stage × season only"))),
                   n = c(ov$n_infection_stage - ov$n_both, ov$n_both, ov$n_three_way_only))
p3a <- ggplot(bars, aes(n, set)) + geom_col(fill = "#526369", width = .6) +
  geom_text(aes(label = comma(n)), hjust = -.15, size = 2.6, family = "Times New Roman") +
  scale_x_continuous(expand = expansion(c(0, .32)), labels = comma, breaks = c(0, 3000, 6000)) +
  labs(title = "Host genes at FDR 0.05", x = "Genes", y = NULL) + theme_fig() + theme(panel.grid.major.y = element_blank())
dr <- fread(file.path(D, "joint_test_permutation_draws.csv"))
obs <- fread(file.path(D, "joint_test_permutation_null.csv"))
TEST <- c(infection_stage = "Infection × stage", infection_stage_season = "Infection × stage × season", only_season_specific = "Three-way only")
long <- melt(dr, id.vars = "draw", measure.vars = names(TEST), variable.name = "test", value.name = "n")
long[, test := factor(TEST[as.character(test)], rev(TEST))]
obs[, test2 := factor(TEST[test], rev(TEST))]
p3b <- ggplot(long, aes(n + 1, test)) +
  geom_point(position = position_jitter(height = .16, width = 0, seed = 2), size = .7, alpha = .45, colour = "#526369", shape = 16) +
  geom_point(data = obs, aes(observed + 1, test2), colour = "#8C2846", shape = 18, size = 3) +
  geom_text(data = obs, aes(observed + 1, test2, label = comma(observed)), vjust = -1, size = 2.5, family = "Times New Roman", colour = "black") +
  scale_x_log10(breaks = c(1, 11, 101, 1001, 10001), labels = c("0", "10", "100", "1,000", "10,000"), limits = c(1, 20000)) +
  labs(title = "Counts under permuted vine labels", subtitle = "Gray: 200 permuted refits. Diamond: observed.", x = "Genes at FDR 0.05 (log scale)", y = NULL) +
  theme_fig() + theme(panel.grid.major.y = element_blank())
at <- fread(file.path(D, "atlas_branch_cell_summary.csv"), colClasses = c(season = "character"))
at <- at[n_genes >= 10]
rank <- at[stage == "Harvest", .(m = max(abs(mean_log2FC_unadjusted))), by = name][order(-m)][1:30]$name
hm <- at[name %in% rank]
ord <- hm[stage == "Harvest", .(m = mean(mean_log2FC_unadjusted)), by = name][order(m)]$name
u <- hm[, .(name, season, col = factor(SHORT[match(stage, STAGES)], SHORT), panel = "Unadjusted", v = mean_log2FC_unadjusted)]
d2 <- hm[stage == "Harvest", .(name, season, col = factor("Har", SHORT), panel = "Adjusted", v = mean_log2FC_adjusted)]
hd <- rbind(u, d2); hd[, name := factor(name, ord)]
hd[, panel := factor(paste(panel, season), c("Unadjusted 2021", "Unadjusted 2022", "Adjusted 2021", "Adjusted 2022"))]
LIM <- 1.5
p3c <- ggplot(hd, aes(col, name, fill = pmax(pmin(v, LIM), -LIM))) + geom_tile(colour = "white", linewidth = .25) +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x", labeller = labeller(panel = function(x) sub(" 20", "\n20", x))) +
  scale_fill_gradient2(low = "#087F8C", mid = "#F4F5F5", high = "#8C2846", limits = c(-LIM, LIM), labels = MINUS, name = "Mean log2FC\n(clipped at ±1.5)") +
  labs(title = "Thirty largest branches at harvest (10 or more genes)", subtitle = "Pre, Ver, Post, Har: stage. Teal lower in infected berries, burgundy higher.", x = NULL, y = NULL) +
  theme_fig() + theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 6.2), legend.position = "right")
save_fig(((p3a | p3b) + plot_layout(widths = c(1, 1.15))) / p3c + plot_layout(heights = c(1, 3.1)) + plot_annotation(tag_levels = "A"), "Figure3", 7, 8.6)

## Figure 4
ex <- fread(file.path(D, "external_stage_concordance.csv"))
ex[, stage := factor(stage, STAGES)]
p4a <- ggplot(ex, aes(stage, spearman_rho, group = 1)) + geom_hline(yintercept = 0, colour = MID, linewidth = .35) +
  geom_line(colour = "#087F8C", linewidth = .6) + geom_point(colour = "#087F8C", size = 2) +
  geom_text(aes(label = gsub("-", "\u2212", sprintf("ρ %.3f\n%.1f%% same sign", spearman_rho, sign_agreement_pct)),
                vjust = ifelse(spearman_rho < 0, 1.35, -.45)), size = 2.5, lineheight = .9, family = "Times New Roman") +
  scale_x_discrete(labels = c("Pre-veraison", "Veraison", "Post-veraison", "Harvest")) + scale_y_continuous(labels = MINUS) + coord_cartesian(ylim = c(-.2, .6)) +
  labs(title = "Merlot against Zinfandel infection effects", subtitle = sprintf("All one-to-one genes (%s at harvest), no significance filter.", comma(ex[stage == "Harvest"]$n_genes)),
       x = NULL, y = "Spearman correlation") + theme_fig()
rows <- data.table(y = 4:1,
  head = c("Shared delay", "Beyond the delay in both seasons", "Beyond the delay mainly in 2022", "Timing and level not separated"),
  body = c("Infected berries sit earlier on the healthy score; larger in 2021 than in 2022",
           "Small heat-shock protein deficit, 27 and 30 of 34 genes",
           "Cytosolic translation down, chloroplast translation up; much of 2021 fits the delay",
           "Photosystem genes earlier than the delay predicts; 11 and 4 of 44 genes remain"))
p4b <- ggplot(rows) + geom_tile(aes(x = .015, y = y), fill = "#526369", width = .006, height = .72) +
  geom_text(aes(x = .035, y = y + .16, label = head), hjust = 0, size = 3, fontface = "bold", family = "Times New Roman") +
  geom_text(aes(x = .035, y = y - .16, label = body), hjust = 0, size = 2.6, family = "Times New Roman", colour = INK) +
  coord_cartesian(xlim = c(0, 1), ylim = c(.5, 4.5), expand = FALSE) +
  labs(title = "Summary of the GRBV response at harvest", subtitle = "Associations in RNA abundance, not measured mechanisms.") + theme_void(base_size = BASE, base_family = "Times New Roman") +
  theme(plot.title = element_text(size = 8.5, face = "bold"), plot.subtitle = element_text(size = 7, colour = INK), plot.tag = element_text(size = 11, face = "bold"), plot.margin = margin(3, 5, 3, 3))
save_fig(p4a / p4b + plot_layout(heights = c(1.1, 1)) + plot_annotation(tag_levels = "A"), "Figure4", 7, 5.6)
cat("T13 COMPLETE\n")
