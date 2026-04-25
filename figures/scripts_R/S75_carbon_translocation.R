# S75 — Carbon translocation / phloem-loading module.  Tests the
# Martinez-Luscher 2018 hexose-translocation hypothesis using a curated 77-gene
# sugar-transport panel (SWEET / SUT / STP / CWINV / VIN / SuSy / SPS /
# hexokinase / trehalose / phloem markers).

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(stringr)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

RESULTS <- "/Users/prempratapsingh/Data_Analysis/farm_pipeline/results/15_carbon_translocation"
OUT_DIR <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S75_carbon_translocation")

# ---- panel A: subclass composition ------------------------------------------
genes <- read_csv(file.path(RESULTS, "sugar_module_genes.csv"), show_col_types = FALSE)
sub_n <- genes %>% count(subclass, name = "n_genes") %>%
  arrange(desc(n_genes)) %>%
  mutate(subclass = factor(subclass, levels = subclass))

pA <- ggplot(sub_n, aes(x = n_genes, y = subclass, fill = subclass)) +
  geom_col(show.legend = FALSE, width = 0.7) +
  geom_text(aes(label = n_genes), hjust = -0.2, size = 3.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "(A) Sugar-transport panel composition",
       subtitle = "77 genes from GFF regex curation",
       x = "# genes", y = NULL) +
  theme_grbv(base_size = 9)

# ---- panel B: module eigengene trajectory -----------------------------------
me <- read_csv(file.path(RESULTS, "sugar_module_eigengene.csv"), show_col_types = FALSE) %>%
  mutate(Maturity_stage = factor(Maturity_stage,
                                 levels = c("Pre_veraison","Veraison","Post_veraison","Harvest"),
                                 labels = c("Pre-ver","Ver","Post-ver","Harvest")),
         Sampling_Season = factor(Sampling_Season))

agg <- me %>%
  group_by(Sampling_Season, Maturity_stage, Infection_Group) %>%
  summarise(mean_me = mean(ME_sugar, na.rm = TRUE),
            se = sd(ME_sugar, na.rm = TRUE) / sqrt(sum(!is.na(ME_sugar))),
            n = n(), .groups = "drop")

pB <- ggplot(agg, aes(x = Maturity_stage, y = mean_me,
                       color = Infection_Group, group = Infection_Group,
                       linetype = Sampling_Season)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
  geom_errorbar(aes(ymin = mean_me - se, ymax = mean_me + se),
                width = 0.15, linewidth = 0.6) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
  scale_color_manual(values = pal_infection, name = NULL) +
  scale_linetype_manual(values = c("2021" = "solid", "2022" = "dashed"), name = NULL) +
  facet_wrap(~ Sampling_Season, ncol = 2) +
  labs(title = "(B) Sugar-transport module eigengene (PC1) across stage x season",
       subtitle = "PC1 = 37.3% of variance on 33 expressed module genes",
       x = NULL, y = "ME_sugar") +
  theme_grbv(base_size = 9) +
  theme(legend.position = "top")

# ---- panel C: per-gene Brix correlation by subclass -------------------------
pg <- read_csv(file.path(RESULTS, "sugar_brix_per_gene.csv"), show_col_types = FALSE) %>%
  filter(!is.na(rho_pooled)) %>%
  mutate(subclass = factor(subclass))

pC <- ggplot(pg, aes(x = subclass, y = rho_pooled, fill = subclass)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, show.legend = FALSE) +
  geom_jitter(width = 0.18, size = 1.6, alpha = 0.7, show.legend = FALSE) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "(C) Per-gene Spearman with harvest Brix",
       subtitle = "Pooled across season; n=27 RB+ samples",
       x = NULL, y = "rho(gene, Brix)") +
  theme_grbv(base_size = 9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

fig <- (pA + pC) / pB + plot_layout(heights = c(0.9, 1)) +
  plot_annotation(title = "S75 — Carbon translocation module",
                  theme = theme(plot.title = element_text(size = 11, face = "bold")))

ggsave(paste0(OUT_STEM, ".png"), fig,
       width = GRBV$W_DOUBLE_MM, height = 200, units = "mm", dpi = 600)
ggsave(paste0(OUT_STEM, ".svg"), fig,
       width = GRBV$W_DOUBLE_MM, height = 200, units = "mm")
cat("[S75] saved\n")
