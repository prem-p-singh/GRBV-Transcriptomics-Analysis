# Main Composite 6 — ML biomarker classification, reframed.
#
# Headline:
#   - All-vars random 5-fold: AUC=1.000 (5 methods) -> "detect virus by viral reads"
#   - Host-only random 5-fold: AUC 0.96-1.00
#   - Host-only LOSO season: AUC 0.76-0.88 (the real generalisation test)
#
# Panels:
#   A  AUC comparison: panel x method x CV scheme (grouped bar)
#   B  AUC drop from random-CV to LOSO season per method (host-only)
#   C  Stability selection: top host genes by LASSO selection frequency
#   D  Feature pool composition (viral vs host)
#
# Inputs:
#   /Users/prempratapsingh/Data_Analysis/advanced_ml_comparison_FIX04.csv
#   /Users/prempratapsingh/Data_Analysis/stability_selection_FIX04.csv
#   /Users/prempratapsingh/Data_Analysis/farm_pipeline/results/fix11_annotations/gene_annotations.csv

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

DATA_ROOT <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR   <- "/Users/prempratapsingh/Data_Analysis/figures/main"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM  <- file.path(OUT_DIR, "Main06_biomarker")

auc <- read_csv(file.path(DATA_ROOT, "advanced_ml_comparison_FIX04.csv"),
                show_col_types = FALSE)
stab <- read_csv(file.path(DATA_ROOT, "stability_selection_FIX04.csv"),
                 show_col_types = FALSE)
ann  <- fread(file.path(DATA_ROOT, "farm_pipeline/results/fix11_annotations/gene_annotations.csv")) %>% as.data.frame()
name_map <- setNames(ann$product, ann$gene_id)

# -----------------------------------------------------------------
# Panel A — AUC comparison
# -----------------------------------------------------------------
auc_long <- auc %>%
  select(panel, method,
         random_5fold_mean_auc, random_5fold_std,
         loco_season_mean_auc,  loco_season_std) %>%
  pivot_longer(cols = c(random_5fold_mean_auc, loco_season_mean_auc),
               names_to = "cv", values_to = "auc") %>%
  mutate(cv = recode(cv,
                     "random_5fold_mean_auc" = "Random 5-fold",
                     "loco_season_mean_auc"  = "Leave-one-season-out"),
         cv = factor(cv, levels = c("Random 5-fold","Leave-one-season-out")),
         panel = recode(panel,
                        "all_vars"  = "All vars (host + viral)",
                        "host_only" = "Host only"))
auc_sd <- auc %>%
  pivot_longer(cols = c(random_5fold_std, loco_season_std),
               names_to = "cv", values_to = "sd") %>%
  mutate(cv = recode(cv,
                     "random_5fold_std" = "Random 5-fold",
                     "loco_season_std"  = "Leave-one-season-out"),
         cv = factor(cv, levels = c("Random 5-fold","Leave-one-season-out")),
         panel = recode(panel,
                        "all_vars"  = "All vars (host + viral)",
                        "host_only" = "Host only"))
auc_long <- auc_long %>% left_join(auc_sd %>% select(panel, method, cv, sd),
                                   by = c("panel","method","cv"))
auc_long$method <- factor(auc_long$method,
                          levels = c("LR","RF","LASSO","SVM-RBF","GB"))

panel_color <- c(
  "All vars (host + viral)" = unname(GRBV_COLORS["RB_pos"]),
  "Host only"               = unname(GRBV_COLORS["RB_neg"])
)

pA <- ggplot(auc_long, aes(x = method, y = auc, fill = panel)) +
  geom_col(position = position_dodge(width = 0.78),
           width = 0.7, color = "white") +
  geom_errorbar(aes(ymin = auc - sd, ymax = auc + sd),
                position = position_dodge(width = 0.78),
                width = 0.20, linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", auc)),
            position = position_dodge(width = 0.78),
            family = "Arial", size = 2.4, vjust = -0.6) +
  facet_wrap(~ cv, ncol = 2) +
  scale_fill_manual(values = panel_color, name = NULL) +
  scale_y_continuous(limits = c(0, 1.10), breaks = seq(0, 1, 0.2),
                     expand = c(0, 0)) +
  geom_hline(yintercept = 0.5, color = "grey60", linetype = "dashed",
             linewidth = 0.3) +
  labs(x = NULL, y = "AUC", tag = "A",
       title = "Classification AUC by panel and CV scheme") +
  theme_grbv() +
  theme(plot.title  = element_text(size = 10, face = "bold"),
        legend.position = "top",
        legend.text  = element_text(size = 7.5),
        strip.text   = element_text(face = "bold", size = 9),
        axis.text.x  = element_text(size = 7.5))

# -----------------------------------------------------------------
# Panel B — AUC drop from random-CV to LOSO (host-only only)
# -----------------------------------------------------------------
drop_df <- auc %>% filter(panel == "host_only") %>%
  transmute(method,
            random = random_5fold_mean_auc,
            loso   = loco_season_mean_auc,
            drop   = random - loso) %>%
  pivot_longer(c(random, loso), names_to = "cv", values_to = "auc") %>%
  mutate(cv = recode(cv,
                     "random" = "Random 5-fold",
                     "loso"   = "LOSO season"),
         method = factor(method, levels = c("LR","RF","LASSO","SVM-RBF","GB")))

drop_label <- drop_df %>% filter(cv == "Random 5-fold") %>%
  left_join(drop_df %>% filter(cv == "LOSO season") %>% select(method, loso = auc),
            by = "method") %>%
  mutate(drop = auc - loso)

pB <- ggplot(drop_df, aes(x = method, y = auc, group = method)) +
  geom_segment(data = drop_label,
               aes(x = method, xend = method,
                   y = auc, yend = loso),
               color = "grey60", linewidth = 0.8) +
  geom_point(aes(color = cv, shape = cv), size = 3.0) +
  geom_text(data = drop_label,
            aes(x = method, y = loso - 0.03,
                label = sprintf("\u0394 = -%.2f", drop)),
            family = "Arial", size = 2.6, color = "grey25",
            hjust = 0.5) +
  scale_color_manual(values = c("Random 5-fold" = unname(GRBV_COLORS["neutral"]),
                                "LOSO season"   = unname(GRBV_COLORS["RB_pos"])),
                     name = NULL) +
  scale_shape_manual(values = c("Random 5-fold" = 16, "LOSO season" = 17),
                     name = NULL) +
  scale_y_continuous(limits = c(0.55, 1.05),
                     breaks = seq(0.6, 1.0, 0.1)) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  labs(x = NULL, y = "AUC (host-only)", tag = "B",
       title = "Generalisation drop:  random CV \u2192 LOSO season") +
  theme_grbv() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = "top",
        legend.text = element_text(size = 7.5),
        axis.text.x = element_text(size = 8))

# -----------------------------------------------------------------
# Panel C — stability selection top host genes
# -----------------------------------------------------------------
stab_top <- stab %>%
  filter(panel == "host_only") %>%
  arrange(desc(lasso_selection_freq)) %>%
  head(20) %>%
  mutate(product = ifelse(gene %in% names(name_map),
                           substr(name_map[gene], 1, 30), ""),
         label = paste0(gene, " \u2014 ", product),
         label = make.unique(label, sep = " #"),
         label = factor(label, levels = rev(label)))

pC <- ggplot(stab_top, aes(x = lasso_selection_freq, y = label)) +
  geom_col(fill = unname(GRBV_COLORS["RB_neg"]), alpha = 0.9, width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", lasso_selection_freq)),
            hjust = -0.18, family = "Arial", size = 2.4, color = "grey20") +
  scale_x_continuous(limits = c(0, 1.05),
                     breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.05))) +
  geom_vline(xintercept = 0.6, color = "grey50",
             linetype = "dashed", linewidth = 0.3) +
  labs(x = "LASSO selection frequency (n=100 subsamples)",
       y = NULL, tag = "C",
       title = "Stability-selected host biomarker candidates") +
  theme_grbv() +
  theme(plot.title  = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 7))

# -----------------------------------------------------------------
# Panel D — feature pool composition
# -----------------------------------------------------------------
n_all  <- unique(auc$n_features[auc$panel == "all_vars"])
n_host <- unique(auc$n_features[auc$panel == "host_only"])
n_viral <- n_all - n_host
pool_df <- data.frame(
  type = c("Viral (N761_*)","Host"),
  n    = c(n_viral, n_host)
) %>% mutate(type = factor(type, levels = c("Viral (N761_*)","Host")),
             pct = n / sum(n))

pD <- ggplot(pool_df, aes(x = "", y = n, fill = type)) +
  geom_col(width = 1, color = "white", linewidth = 0.6) +
  coord_polar(theta = "y") +
  geom_text(aes(label = sprintf("%s\n(n = %d, %.1f%%)",
                                type, n, 100 * pct)),
            position = position_stack(vjust = 0.5),
            family = "Arial", size = 2.7, fontface = "bold",
            color = ifelse(pool_df$type == "Viral (N761_*)", "white", "grey15")) +
  scale_fill_manual(values = c("Viral (N761_*)" = unname(GRBV_COLORS["RB_pos"]),
                               "Host"           = unname(GRBV_COLORS["RB_neg"])),
                    guide = "none") +
  labs(tag = "D", title = "Feature pool composition",
       caption = sprintf("Top %d most variable genes;  n = %d viral, %d host",
                         n_all, n_viral, n_host)) +
  theme_void(base_family = "Arial") +
  theme(plot.title   = element_text(size = 10, face = "bold", hjust = 0.5,
                                    family = "Arial"),
        plot.caption = element_text(size = 7, color = "grey35", hjust = 0.5,
                                    margin = margin(t = 4), family = "Arial"),
        plot.tag     = element_text(face = "bold", size = 13, family = "Arial"),
        plot.tag.position = c(0.02, 0.98))

# -----------------------------------------------------------------
# Composite layout
# -----------------------------------------------------------------
layout <- pA / (pB | pD) / pC +
  plot_layout(heights = c(1.0, 0.95, 1.4),
              widths  = c(1, 1)) +
  plot_annotation(
    title = "Reframed ML biomarker analysis: AUC=1.000 was an artifact of viral reads",
    subtitle = "Random 5-fold AUC stays \u2248 1.000 with all variables, drops to 0.76-0.88 under leave-one-season-out CV with host genes only",
    theme = theme(plot.title = element_text(family = "Arial", face = "bold",
                                            size = 13, hjust = 0.5,
                                            margin = margin(b = 2)),
                  plot.subtitle = element_text(family = "Arial", size = 8.5,
                                               hjust = 0.5, color = "grey25",
                                               margin = margin(b = 6)))
  )

# Widen composite so panel B has room for 5 methods + \u0394 labels and
# panel C can display long "LOC######## \u2014 <product>" gene labels.
save_grbv(layout, OUT_STEM,
          width_mm = 220, height_mm = 280)

cat("[main06] AUC summary:\n"); print(auc_long, n = 20)
