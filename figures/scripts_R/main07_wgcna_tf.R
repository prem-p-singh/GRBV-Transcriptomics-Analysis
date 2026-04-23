# Main Composite 7 — WGCNA modules + transcription-factor regulators.
#
# Panels:
#   A  Module-trait correlation heatmap (13 modules x 4 traits)
#   B  TF family DEG enrichment heatmap across the 8 stage x season conditions
#   C  Hub genes (top kME) for the most infection-correlated module
#   D  TF family mean log2FC trajectory across stages (split by season)
#
# Inputs:
#   /Users/prempratapsingh/Data_Analysis/module_trait_correlations.csv
#   /Users/prempratapsingh/Data_Analysis/coexpression_modules.csv
#   /Users/prempratapsingh/Data_Analysis/hub_genes_per_module.csv
#   /Users/prempratapsingh/Data_Analysis/tf_family_DE_summary.csv

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/.claude/worktrees/gallant-sinoussi-04711f/farm_pipeline/scripts/figures"
source(file.path(HERE, "theme_grbv.R"))

DATA_ROOT <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR   <- "/Users/prempratapsingh/Data_Analysis/figures/main"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM  <- file.path(OUT_DIR, "Main07_wgcna_tf")

mt <- fread(file.path(DATA_ROOT, "module_trait_correlations.csv")) %>% as.data.frame()
colnames(mt)[1] <- "module"
mod <- fread(file.path(DATA_ROOT, "coexpression_modules.csv")) %>% as.data.frame()
hub <- fread(file.path(DATA_ROOT, "hub_genes_per_module.csv")) %>% as.data.frame()
tf  <- fread(file.path(DATA_ROOT, "tf_family_DE_summary.csv")) %>% as.data.frame()

# -----------------------------------------------------------------
# Panel A — module-trait correlation heatmap
# -----------------------------------------------------------------
mt_long <- mt %>%
  pivot_longer(-module, names_to = "trait", values_to = "r") %>%
  mutate(trait = recode(trait,
                        "InfectionYear" = "Years infected"))
# Order modules by abs(r) for Infection
ord <- mt %>% mutate(abs_inf = abs(Infection)) %>%
  arrange(desc(abs_inf)) %>% pull(module)
mt_long$module <- factor(mt_long$module, levels = rev(ord))
mt_long$trait  <- factor(mt_long$trait,
                         levels = c("Stage","Season","Infection","Years infected"))

pA <- ggplot(mt_long, aes(x = trait, y = module, fill = r)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(abs(r) < 0.005, "0",
                               formatC(r, format = "f", digits = 2,
                                       flag = "+"))),
            family = "Arial", size = 2.0,
            color = ifelse(abs(mt_long$r) > 0.5, "white", "grey15")) +
  scale_fill_gradient2(low = unname(GRBV_COLORS["down"]),
                       mid = "white",
                       high = unname(GRBV_COLORS["up"]),
                       midpoint = 0, limits = c(-1, 1),
                       name = "Pearson r") +
  labs(x = NULL, y = NULL, tag = "A",
       title = "Module-trait correlations") +
  theme_grbv() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        legend.key.height = unit(20, "pt"))

# -----------------------------------------------------------------
# Panel B — TF family DEG enrichment heatmap
# -----------------------------------------------------------------
tf <- tf %>%
  separate(condition, into = c("season","stage_a","stage_b"),
           sep = "_", remove = FALSE, fill = "right") %>%
  mutate(stage = ifelse(is.na(stage_b), stage_a,
                        paste(stage_a, stage_b, sep = "_")),
         stage = factor(stage, levels = GRBV_STAGES),
         season = factor(season, levels = c("2021","2022")),
         tf_family = factor(tf_family,
                            levels = c("ERF/AP2","MYB","WRKY","bHLH","NAC","HSF")))

pB <- ggplot(tf, aes(x = stage, y = tf_family,
                     fill = percent_de,
                     label = ifelse(percent_de < 0.1, "",
                                    sprintf("%.0f", percent_de)))) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(color = percent_de > 12),
            family = "Arial", size = 2.0) +
  facet_wrap(~ season, ncol = 2) +
  scale_x_discrete(labels = GRBV_STAGE_LABEL) +
  scale_fill_gradient(low = "#fffaf0", high = unname(GRBV_COLORS["RB_pos"]),
                      name = "% DE") +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"),
                     guide = "none") +
  labs(x = NULL, y = NULL, tag = "B",
       title = "TF family DEG enrichment") +
  theme_grbv() +
  theme(plot.title  = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 7.5),
        axis.text.y = element_text(size = 8),
        strip.text  = element_text(face = "bold", size = 9),
        legend.key.height = unit(15, "pt"))

# -----------------------------------------------------------------
# Panel C — top hub genes for the most infection-correlated module
# -----------------------------------------------------------------
target_mod <- ord[1]  # the |r| top module on Infection
hub_top <- hub %>%
  filter(module == target_mod) %>%
  arrange(desc(kME)) %>% head(10) %>%
  mutate(label = paste0(gene, " \u2014 ", substr(gene_name, 1, 28)),
         label = make.unique(label, sep = " #"),
         label = factor(label, levels = rev(label)))

pC <- ggplot(hub_top, aes(x = kME, y = label, fill = log2FC)) +
  geom_col(alpha = 0.92, width = 0.7, color = "white") +
  geom_text(aes(label = sprintf("kME=%.2f, FC=%+.2f", kME, log2FC)),
            hjust = -0.10, family = "Arial", size = 2.3, color = "grey15") +
  scale_fill_gradient2(low = unname(GRBV_COLORS["down"]),
                       mid = "white",
                       high = unname(GRBV_COLORS["up"]),
                       midpoint = 0, name = "log2FC", limits = c(-3, 3),
                       oob = scales::squish) +
  scale_x_continuous(limits = c(0, 1.05),
                     breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.75))) +
  labs(x = "Module membership (kME)", y = NULL, tag = "C",
       title = sprintf("Hub genes: '%s' module", target_mod)) +
  theme_grbv() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 7),
        legend.key.height = unit(12, "pt"))

# -----------------------------------------------------------------
# Panel D — TF family mean log2FC per stage (line trajectories)
# -----------------------------------------------------------------
tf_lines <- tf %>%
  group_by(season, stage, tf_family) %>%
  summarise(mean_FC = mean(mean_log2FC, na.rm = TRUE), .groups = "drop") %>%
  mutate(stage = factor(stage, levels = GRBV_STAGES))

family_color <- c("ERF/AP2" = "#e74c3c", "MYB" = "#3498db",
                  "WRKY"   = "#9b59b6", "bHLH" = "#2ecc71",
                  "NAC"    = "#f39c12", "HSF"  = "#16a085")

pD <- ggplot(tf_lines, aes(x = stage, y = mean_FC,
                           color = tf_family, group = tf_family)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.0) +
  facet_wrap(~ season, ncol = 2) +
  scale_color_manual(values = family_color, name = "TF family") +
  scale_x_discrete(labels = GRBV_STAGE_LABEL) +
  labs(x = NULL,
       y = "Mean log2FC (RB+ vs RB\u2212)",
       tag = "D",
       title = "TF family FC trajectories") +
  theme_grbv() +
  theme(plot.title  = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 7.5),
        strip.text  = element_text(face = "bold", size = 9),
        legend.position = "right",
        legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 7.5, face = "bold"),
        legend.key.height = unit(7, "pt"))

# -----------------------------------------------------------------
# Composite layout
# -----------------------------------------------------------------
layout <- (pA | pB) / (pC | pD) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Co-expression modules and transcription-factor regulators",
    subtitle = "13 WGCNA modules across 4 traits;  6 TF families (267 genes) profiled across 8 stage \u00d7 season conditions",
    theme = theme(plot.title = element_text(family = "Arial", face = "bold",
                                            size = 13, hjust = 0.5,
                                            margin = margin(b = 2)),
                  plot.subtitle = element_text(family = "Arial", size = 8.5,
                                               hjust = 0.5, color = "grey25",
                                               margin = margin(b = 6)))
  )

# Widen composite so panels B + D (each with 2 stage-labelled facets) and
# panel C (long "kME=..., FC=..." annotations on the right) do not clip.
save_grbv(layout, OUT_STEM,
          width_mm = 230, height_mm = 230)

cat(sprintf("[main07] target module for hub panel: '%s'\n", target_mod))
cat(sprintf("[main07] modules with |r_infection|>0.10: %d\n",
            sum(abs(mt$Infection) > 0.10)))
