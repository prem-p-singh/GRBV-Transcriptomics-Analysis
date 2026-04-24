# S72 — Per-gene variance-decomposition ternary plot.
#
# Conceptual replica of Singh et al. 2025 bioRxiv Figure 5 (B. cinerea host x
# isolate variance partition), adapted to the GRBV/Merlot single-species
# design.  Axis translation:
#
#   Singh 2025 ternary       Our ternary (same math, different axes)
#   -----------------------  -----------------------------------------
#   Host species (H\u00b2_host)    -> Maturity stage (eta^2_Stg)
#   Pathogen isolate (H\u00b2_iso) -> Infection + (Inf x Stg) + (Inf x Sea)
#                                 + (Inf x Sea x Stg)
#   Host x Isolate (H\u00b2_HxI)   -> Season + (Sea x Stg) non-infection env
#
# That is: axis 1 is developmental, axis 2 is viral, axis 3 is environmental.
# Each gene's three eta^2 contributions are renormalised so their coordinates
# sum to 1 on the simplex (the rest of the variance is residual/noise).
#
# Output: figures/supp/S72_variance_decomposition_ternary.png + .svg

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table); library(stringr)
  library(ggtern)
})

HERE     <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ROOT     <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR  <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S72_variance_decomposition_ternary")

# ============================================================================
# LOAD
# ============================================================================
anova <- fread(file.path(ROOT, "multifactor_anova_results.csv")) %>%
  as.data.frame()

# Aggregate eta^2 onto three semantic buckets
anova$Developmental <- anova$eta2_Stg
anova$Viral         <- anova$eta2_Inf +
                       anova$`eta2_Inf:Stg` +
                       anova$`eta2_Inf:Sea` +
                       anova$`eta2_Inf:Sea:Stg`
anova$Environmental <- anova$eta2_Sea +
                       anova$`eta2_Sea:Stg`

# Renormalise each gene to the simplex: rescale so Dev+Vir+Env = 1.
# Keep the original total as "total_explained" for diagnostics.
anova$total_explained <- with(anova, Developmental + Viral + Environmental)
# Drop genes where no component explains > 0 (residual-dominated)
anova <- anova %>% filter(total_explained > 0.02)

anova <- anova %>%
  mutate(Developmental = Developmental / total_explained,
         Viral         = Viral         / total_explained,
         Environmental = Environmental / total_explained)

cat(sprintf("[S72] retained %d of 14,604 genes with total explained >2%%\n",
            nrow(anova)))

# ============================================================================
# Overlay categories
# ============================================================================
heat_de <- fread(file.path(ROOT, "heat_stress_DEGs_all_conditions.csv")) %>%
  as.data.frame()
shsp_genes <- unique(heat_de$Gene[grepl("sHSP", heat_de$Category)])
hsf_genes  <- unique(heat_de$Gene[grepl("HSF",  heat_de$Category)])

tf_fam <- fread(file.path(ROOT, "tf_family_gene_details.csv")) %>% as.data.frame()
wrky_genes <- unique(tf_fam$gene[tf_fam$tf_family == "WRKY"])

sig <- read_csv(file.path(ROOT, "tier1_grbv_signature.csv"),
                show_col_types = FALSE) %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FC) >= 1)
grbv_sig_genes <- sig$LOC_id

# Tag each gene with its category (one category only; priority order)
anova$category <- "other"
anova$category[anova$gene %in% wrky_genes]    <- "WRKY TF"
anova$category[anova$gene %in% hsf_genes]     <- "HSF"
anova$category[anova$gene %in% shsp_genes]    <- "sHSP"
anova$category[anova$gene %in% grbv_sig_genes]<- "GRBV signature"
anova$category <- factor(anova$category,
                         levels = c("other","WRKY TF","HSF","sHSP","GRBV signature"))

cat_counts <- anova %>% count(category)
print(cat_counts)

cat_pal <- c("other"          = alpha("grey75", 0.10),
             "WRKY TF"        = unname(GRBV_COLORS["season_2022"]),
             "HSF"            = unname(GRBV_COLORS["warm"]),
             "sHSP"           = unname(GRBV_COLORS["RB_pos"]),
             "GRBV signature" = unname(GRBV_COLORS["med_titer"]))

# ============================================================================
# Panel A — ggtern density of all genes
# ============================================================================
# Shared theme block applied to BOTH ternary panels so they have identical
# proportions and border weight.  Uses ggtern's native axis.line elements
# (no manual overlay), with linewidth syntax (ggplot >= 3.4) which ggtern
# 4.0 understands.  coord_tern() is used with defaults (no custom limits).
# Font sizes here are ~2.5x panel C's native sizes because the intermediate
# 9.5"x8" raster is scaled down by that factor when embedded in the composite.
# On the final page, these translate to ~11 / 10 / 10 / 8 pt -- matching
# panel C's 10 / 8 / 9 / 7.5 pt footprint.
tern_theme <- theme_bw(base_family = "Arial", base_size = 26) +
  theme_showarrows() +
  theme(plot.title    = element_text(size = 28, face = "bold", family = "Arial"),
        plot.subtitle = element_text(size = 22, colour = "grey30", family = "Arial"),
        tern.axis.title.T = element_text(size = 24, face = "bold"),
        tern.axis.title.L = element_text(size = 24, face = "bold"),
        tern.axis.title.R = element_text(size = 24, face = "bold"),
        tern.axis.text  = element_text(size = 18),
        tern.axis.line.T = element_line(colour = "black", linewidth = 1.6),
        tern.axis.line.L = element_line(colour = "black", linewidth = 1.6),
        tern.axis.line.R = element_line(colour = "black", linewidth = 1.6),
        tern.axis.ticks = element_line(colour = "black", linewidth = 0.6),
        tern.panel.grid.minor = element_blank(),
        tern.panel.grid.major = element_line(colour = "grey88", linewidth = 0.5),
        tern.panel.background = element_blank(),
        tern.plot.background  = element_blank(),
        legend.position   = "bottom",
        legend.key.height = unit(24, "pt"),
        legend.key.width  = unit(26, "pt"),
        legend.title = element_text(size = 22),
        legend.text  = element_text(size = 20),
        plot.margin  = margin(20, 40, 20, 40))

# Panel A: ternary point scatter of ALL 14,604 genes.  Same presentation
# pattern as panel B (geom_point on ggtern) but without the category
# overlay — natural density emerges from overplotting with low alpha.
pA <- ggtern(data = anova,
             aes(x = Developmental, y = Viral, z = Environmental)) +
  geom_point(colour = unname(GRBV_COLORS["stage_pre"]),
             alpha = 0.12,
             size = 0.45) +
  Llab("Dev") + Tlab("Viral") + Rlab("Env") +
  labs(title = "A  Gene-level variance decomposition",
       subtitle = sprintf("%s genes \u00b7 density via overplotting",
                          scales::comma(nrow(anova)))) +
  tern_theme

# ============================================================================
# Panel B — categories overlaid on ternary
# ============================================================================
anova_plot <- anova %>% arrange(category)  # draw "other" first

# Panel B: categories overlaid — same tern_theme so borders match A exactly.
pB <- ggtern(data = anova_plot,
             aes(x = Developmental, y = Viral, z = Environmental,
                 colour = category, size = category, alpha = category)) +
  geom_point() +
  scale_colour_manual(values = cat_pal, name = NULL,
                      breaks = c("GRBV signature","sHSP","HSF","WRKY TF","other")) +
  scale_size_manual(values = c("other" = 0.45,
                               "WRKY TF" = 1.3,
                               "HSF" = 1.3,
                               "sHSP" = 1.4,
                               "GRBV signature" = 1.4),
                    guide = "none") +
  scale_alpha_manual(values = c("other" = 0.15,
                                "WRKY TF" = 0.85,
                                "HSF" = 0.90,
                                "sHSP" = 0.90,
                                "GRBV signature" = 0.90),
                     guide = "none") +
  Llab("Dev") + Tlab("Viral") + Rlab("Env") +
  labs(title = "B  Biological categories in variance space",
       subtitle = "GRBV signature \u2192 Viral; sHSPs + HSFs \u2192 Dev-Env") +
  tern_theme

# ============================================================================
# Panel C — per-category eta^2 box plots (quantitative view)
# ============================================================================
box_df <- anova %>%
  filter(category != "other") %>%
  select(category, Developmental, Viral, Environmental) %>%
  pivot_longer(c(Developmental, Viral, Environmental),
               names_to = "axis", values_to = "value") %>%
  mutate(axis = factor(axis,
                       levels = c("Developmental","Viral","Environmental")))

axis_pal <- c("Developmental" = unname(GRBV_COLORS["stage_pre"]),
              "Viral"         = unname(GRBV_COLORS["RB_pos"]),
              "Environmental" = unname(GRBV_COLORS["season_2021"]))

pC <- ggplot(box_df, aes(x = axis, y = value, fill = axis)) +
  geom_boxplot(outlier.size = 0.4, width = 0.65, colour = "grey25",
               linewidth = 0.3) +
  facet_wrap(~ category, ncol = 4) +
  scale_fill_manual(values = axis_pal, guide = "none") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = NULL, y = "Normalised \u03b7\u00b2 share (sums to 1 per gene)",
       title = "C  Mean variance share per category",
       subtitle = "GRBV signature: predominantly viral; sHSPs: split viral + developmental; WRKY: environmental tail") +
  theme_grbv(base_size = 9) +
  theme(plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        strip.text    = element_text(face = "bold", size = 9),
        axis.text.x   = element_text(angle = 25, hjust = 1, size = 7.5),
        axis.text.y   = element_text(size = 7.5),
        panel.spacing = unit(6, "pt"))

# ============================================================================
# Compose: ggtern panels need grid assembly (patchwork + ggtern is fiddly,
# so use gridExtra::arrangeGrob for the top row, then stack with plain ggplot)
# ============================================================================
suppressPackageStartupMessages({library(gridExtra); library(grid); library(png)})

# ggtern's tern.axis.line.T/L/R theme elements are dropped by ggplotGrob()
# during arrangeGrob composition, so the triangle border disappears.  Fix:
# render each ternary panel to a temp PNG via direct grid.draw (which keeps
# the axis-line elements intact), then re-import as rasterGrobs.  pC is a
# regular ggplot, composed normally.
render_tern <- function(plot, stem, w_in = 9.5, h_in = 8.0, res = 600) {
  fout <- tempfile(pattern = stem, fileext = ".png")
  ragg::agg_png(filename = fout, width = w_in, height = h_in,
                units = "in", res = res, background = "white")
  grid.draw(ggplotGrob(plot))
  invisible(dev.off())
  rasterGrob(readPNG(fout), interpolate = TRUE)
}
pA_raster <- render_tern(pA, "tern_pA")
pB_raster <- render_tern(pB, "tern_pB")

top_row <- arrangeGrob(pA_raster, pB_raster, ncol = 2, widths = c(1, 1))
pC_grob <- ggplotGrob(pC)

final <- arrangeGrob(top_row, pC_grob, ncol = 1,
                     heights = c(1.3, 0.7))

# Add a title block
title_grob <- textGrob(
  "Fig. S72 \u2014 Gene-level variance decomposition (developmental / viral / environmental)",
  gp = gpar(fontface = "bold", fontsize = 11, fontfamily = "Arial"))
subtitle_grob <- textGrob(
  "Conceptual replica of Singh et al. 2025 bioRxiv Fig 5 ternary, adapted to the GRBV/Merlot single-species design.",
  gp = gpar(fontsize = 8.5, col = "grey25", fontfamily = "Arial"))

full <- arrangeGrob(title_grob, subtitle_grob, final,
                    ncol = 1, heights = c(0.04, 0.03, 0.93))

# ggsave cannot scale raw grobs correctly — use an explicit device
w_in <- GRBV$W_DOUBLE_MM / 25.4
h_in <- 235 / 25.4

ragg::agg_png(filename = paste0(OUT_STEM, ".png"),
              width = w_in, height = h_in, units = "in",
              res = 600, background = "white")
grid::grid.draw(full)
invisible(dev.off())

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(filename = paste0(OUT_STEM, ".svg"),
                   width = w_in, height = h_in, bg = "white")
  grid::grid.draw(full)
  invisible(dev.off())
}

cat(sprintf("[S72] wrote %s.png / .svg\n", OUT_STEM))

# ============================================================================
# Summary statistics for Discussion
# ============================================================================
summary_table <- anova %>% group_by(category) %>%
  summarise(n = n(),
            dev_mean = mean(Developmental),
            vir_mean = mean(Viral),
            env_mean = mean(Environmental),
            .groups = "drop")
print(summary_table)
