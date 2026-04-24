# S70 — Literature-curated GRBV-related pathway enrichment across conditions.
#
# Filters the existing GSEA panel (375 gene sets x 8 season x stage conditions)
# down to the ~30 pathways reported in the GRBV / Vitis-virus literature to
# be affected by grapevine Red Blotch Virus infection, then shows how NES
# changes across the season x stage grid.
#
# Literature basis for the curated category list (coarse):
#   Sudarshana et al. 2015  Annu Rev Phytopathol  : sugar / ripening blockade
#   Blanco-Ulate et al. 2017 J Exp Bot 68:1225    : anthocyanin, phenylpropanoid,
#                                                   hormone signaling, sugar
#   Yepes et al. 2018 Planta 248:395              : photosynthesis, carbohydrate,
#                                                   hormone signaling
#   Chitarra et al. 2018 BMC Plant Biol 18:41     : secondary metabolism,
#                                                   flavonoid
#   Song et al. 2022 Viruses 14:1831              : heat shock, ribosome,
#                                                   protein folding
#   Agudelo-Romero et al. 2015 J Exp Bot 66:1769  : pathogen response,
#                                                   hormone signaling
#   This study (Tier 1 meta-analysis)             : GRBV signature overlap
#
# Output: figures/supp/S70_literature_curated_pathways.png + .svg

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table); library(stringr)
  library(forcats)
})

HERE     <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ROOT     <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR  <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S70_literature_curated_pathways")

STAGES    <- c("Pre_veraison","Veraison","Post_veraison","Harvest")
STAGE_LAB <- c(Pre_veraison = "Pre-v\u00e9r", Veraison = "V\u00e9r",
               Post_veraison = "Post-v\u00e9r", Harvest = "Harvest")
COND_ORDER <- c("2021_Pre_veraison","2021_Veraison",
                "2021_Post_veraison","2021_Harvest",
                "2022_Pre_veraison","2022_Veraison",
                "2022_Post_veraison","2022_Harvest")

# ============================================================================
# LITERATURE-CURATED CATEGORY TABLE
# Regex patterns map our GSEA gene_set names to categories that the GRBV
# and Vitis-virus literature has flagged as relevant.
# ============================================================================
lit_table <- tribble(
  ~category,                   ~pattern,                                             ~direction_lit, ~ref,
  "Sugar / starch",            "[Ss]tarch|[Ss]ucrose|[Gg]lycolysis|[Gg]luconeogen|[Ff]ructose|[Hh]exose|[Nn]ucleotide sugar|[Pp]entose phosphate", "mixed",    "Sudarshana 15; Blanco-Ulate 17",
  "Phenylpropanoid / flavonoid", "[Pp]henylpropan|[Ff]lavonoid|[Ss]econdary metabolite|[Tt]erpenoid|[Cc]arotenoid|[Ii]soquinoline",               "down",      "Blanco-Ulate 17; Chitarra 18",
  "Photosynthesis / chloroplast", "[Pp]hotosynth|[Pp]hotosystem|light harvesting|[Cc]hloroplast|[Tt]hylakoid|[Cc]alvin|[Pp]lastid|[Cc]arbon fixation|[Cc]hlorophyll", "down", "Yepes 18; Sudarshana 15",
  "Heat shock / chaperone",    "heat shock|[Cc]haperone|HSP|[Uu]nfolded protein|[Pp]rotein folding|[Pp]rotein complex oligomer|[Pp]rotein processing in endoplasmic", "down", "This study; Song 22",
  "Plant hormone signaling",   "[Pp]lant hormone|[Aa]bscisic|[Ee]thylene|[Aa]uxin|[Jj]asmon|[Ss]alicyl|[Hh]ormone-mediated",                   "mixed",     "Blanco-Ulate 17; Agudelo 15",
  "Oxidative / ROS",           "[Gg]lutathione|[Oo]xidative|hydrogen peroxide|[Pp]eroxid|[Ss]uperoxide|[Aa]scorbat|[Rr]eactive oxygen|[Oo]xidored", "up",      "Agudelo 15; Blanco-Ulate 17",
  "Defense / immunity",        "[Dd]efense|[Dd]efence|[Ii]mmune|[Pp]athogenesis|[Rr]esponse to bacter|[Rr]esponse to fungus",                          "mixed",     "Fung 08; Polesani 10; Song 22",
  "Ribosome / translation",    "^KEGG:[Rr]ibosom|[Tt]ranslation initiation|[Tt]ranslation|ribosomal subunit|rRNA|small subunit|large subunit",  "down",     "This study; Song 22",
  "Cell wall / pectin",        "[Cc]ell wall|[Pp]ectin|[Cc]ellulose|xyloglucan|glucan|hydrolyz|apoplast",                                      "mixed",     "Blanco-Ulate 17; Agudelo 15",
  "Amino acid metabolism",     "[Aa]mino acid|[Aa]rginin|[Ll]ysin|[Pp]roline|[Gg]lutam|[Ss]erine|[Ll]eucin|[Vv]aline|[Aa]lanin|[Aa]spartat|[Ll]ysosome|[Pp]roteasome", "mixed", "Yepes 18; Blanco-Ulate 17"
)

# ============================================================================
# LOAD + MATCH
# ============================================================================
gs <- read_csv(file.path(ROOT, "gsea_results_all_conditions.csv"),
               show_col_types = FALSE) %>%
  mutate(Season = str_extract(condition, "^\\d{4}"),
         Stage  = str_remove(condition, "^\\d{4}_"),
         Stage  = factor(Stage, levels = STAGES),
         condition = factor(condition, levels = COND_ORDER))

# Assign each gene_set to the FIRST matching literature category (one-to-one)
assign_cat <- function(name) {
  for (i in seq_len(nrow(lit_table))) {
    if (str_detect(name, lit_table$pattern[i])) return(lit_table$category[i])
  }
  return(NA_character_)
}
gs$category <- vapply(gs$gene_set, assign_cat, character(1))

# Restrict to literature-matched sets
gs_lit <- gs %>% filter(!is.na(category)) %>%
  mutate(category = factor(category, levels = lit_table$category))

cat(sprintf("[S70] %d of 375 GSEA pathways matched literature categories\n",
            length(unique(gs_lit$gene_set))))

# Per-category mean NES across conditions
cat_summary <- gs_lit %>%
  group_by(category, condition, Season, Stage) %>%
  summarise(mean_nes  = mean(nes, na.rm = TRUE),
            n_sets    = n_distinct(gene_set),
            any_sig   = any(fdr < 0.05, na.rm = TRUE),
            .groups = "drop")

# Count per category
cat_count <- gs_lit %>% group_by(category) %>%
  summarise(n = n_distinct(gene_set), .groups = "drop")

# ============================================================================
# Panel A — Category-level mean NES heatmap (rows = lit category)
# ============================================================================
cat_summary$cell_label <- sprintf("%+.2f%s",
                                   cat_summary$mean_nes,
                                   ifelse(cat_summary$any_sig, "*", ""))

lim_A <- max(abs(cat_summary$mean_nes), na.rm = TRUE)
lim_A <- max(lim_A, 1.5)

pA <- ggplot(cat_summary,
             aes(x = condition, y = fct_rev(category), fill = mean_nes)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = cell_label),
            family = "Arial", size = 2.4, colour = "grey15") +
  scale_fill_gradient2(low = unname(GRBV_COLORS["down"]),
                       mid = "white",
                       high = unname(GRBV_COLORS["up"]),
                       midpoint = 0, limits = c(-lim_A, lim_A),
                       oob = scales::squish, name = "mean NES") +
  scale_x_discrete(labels = function(x) {
    xs <- str_split_fixed(x, "_", 2)
    paste0(xs[, 1], "\n", STAGE_LAB[xs[, 2]])
  }) +
  labs(x = NULL, y = NULL,
       title = "A  Literature-curated category NES across conditions",
       subtitle = "mean NES across sets in each category; * = at least one set FDR<0.05 in that cell") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        axis.text.x   = element_text(size = 7, lineheight = 1.0),
        axis.text.y   = element_text(size = 8),
        legend.key.width  = unit(8, "pt"),
        legend.key.height = unit(14, "pt"))

# ============================================================================
# Panel B — Individual pathway heatmap (grouped by category)
# ============================================================================
# Pick up to 4 most variable pathways per category
top_per_cat <- gs_lit %>%
  group_by(gene_set, category) %>%
  summarise(max_nes = max(abs(nes), na.rm = TRUE), .groups = "drop") %>%
  group_by(category) %>%
  slice_max(max_nes, n = 4, with_ties = FALSE) %>%
  ungroup()

matB <- gs_lit %>% filter(gene_set %in% top_per_cat$gene_set)

# Short label with category prefix tag
matB <- matB %>%
  mutate(cat_short = case_when(
    category == "Sugar / starch"            ~ "SUG",
    category == "Phenylpropanoid / flavonoid" ~ "FLA",
    category == "Photosynthesis / chloroplast" ~ "PHO",
    category == "Heat shock / chaperone"    ~ "HSP",
    category == "Plant hormone signaling"   ~ "HOR",
    category == "Oxidative / ROS"           ~ "OXR",
    category == "Defense / immunity"        ~ "DEF",
    category == "Ribosome / translation"    ~ "RIB",
    category == "Cell wall / pectin"        ~ "CWL",
    category == "Amino acid metabolism"     ~ "AA ",
    TRUE                                    ~ "???"),
    gs_short = str_remove(gene_set, "^[A-Z]+:"),
    gs_short = ifelse(nchar(gs_short) > 32,
                      paste0(substr(gs_short, 1, 31), "\u2026"),
                      gs_short),
    gs_tag   = paste0("[", cat_short, "] ", gs_short))

# Row order: by category then by pathway name
row_order_B <- matB %>% arrange(category, gene_set) %>%
  distinct(gs_tag) %>% pull(gs_tag) %>% rev()
matB$gs_tag <- factor(matB$gs_tag, levels = row_order_B)

pB <- ggplot(matB, aes(x = condition, y = gs_tag, fill = nes)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(fdr < 0.05, "*", "")),
            family = "Arial", size = 2.4, colour = "grey15", vjust = 0.75) +
  scale_fill_gradient2(low = unname(GRBV_COLORS["down"]),
                       mid = "white",
                       high = unname(GRBV_COLORS["up"]),
                       midpoint = 0, limits = c(-3, 3),
                       oob = scales::squish, name = "NES") +
  scale_x_discrete(labels = function(x) {
    xs <- str_split_fixed(x, "_", 2)
    paste0(xs[, 1], "\n", STAGE_LAB[xs[, 2]])
  }) +
  labs(x = NULL, y = NULL,
       title = "B  Top 4 pathways per category, individual NES",
       subtitle = "rows grouped by literature category (tag in brackets); * = FDR<0.05") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        axis.text.x   = element_text(size = 6.5, lineheight = 1.0),
        axis.text.y   = element_text(size = 6.5),
        legend.key.width  = unit(8, "pt"),
        legend.key.height = unit(14, "pt"))

# ============================================================================
# Panel C — Per-category NES trajectory lines (stage x season)
# ============================================================================
# Shorter category aliases so 10 facet strips don't clip
short_cat_map <- c(
  "Sugar / starch"                = "Sugar",
  "Phenylpropanoid / flavonoid"   = "Flavonoid",
  "Photosynthesis / chloroplast"  = "Photosynth.",
  "Heat shock / chaperone"        = "HSP / chap.",
  "Plant hormone signaling"       = "Hormone",
  "Oxidative / ROS"               = "Oxidative",
  "Defense / immunity"            = "Defense",
  "Ribosome / translation"        = "Ribosome",
  "Cell wall / pectin"            = "Cell wall",
  "Amino acid metabolism"         = "AA metab."
)
cat_summary$cat_short <- factor(short_cat_map[as.character(cat_summary$category)],
                                 levels = unname(short_cat_map))

pC <- ggplot(cat_summary, aes(x = Stage, y = mean_nes,
                               colour = Season, group = Season)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = c("2021" = unname(GRBV_COLORS["season_2021"]),
                                 "2022" = unname(GRBV_COLORS["season_2022"])),
                      name = "Season") +
  scale_x_discrete(labels = STAGE_LAB) +
  facet_wrap(~ cat_short, ncol = 5, scales = "free_y") +
  labs(x = NULL, y = "Mean NES",
       title = "C  Pathway trajectory per literature category",
       subtitle = "mean NES across all sets in each category, plotted per season") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        strip.text    = element_text(face = "bold", size = 7.5),
        axis.text.x   = element_text(angle = 25, hjust = 1, size = 7),
        axis.text.y   = element_text(size = 7),
        legend.position = "top",
        panel.spacing.x = unit(6, "pt"),
        panel.spacing.y = unit(4, "pt"))

# ============================================================================
# Compose
# ============================================================================
# Streamlined 2-panel layout: category-level heatmap (A) + trajectories (C).
# Panel B (individual pathway detail) is exported as a companion CSV
# instead of a sub-panel so this figure stays readable.
matB_out <- matB %>% select(category, gene_set, condition, nes, fdr) %>%
  arrange(category, gene_set, condition)
write.csv(matB_out,
          "/Users/prempratapsingh/Data_Analysis/S70_top4_per_category_NES.csv",
          row.names = FALSE)

layout <- pA / pC +
  plot_layout(heights = c(1.1, 1.0)) +
  plot_annotation(
    title = "Fig. S70 \u2014 GRBV-related pathway enrichment shifts (literature-curated)",
    subtitle = sprintf("Filtered to %d of 375 GSEA pathways matching 10 categories flagged by the GRBV / Vitis-virus literature (Sudarshana 15, Blanco-Ulate 17, Yepes 18, Chitarra 18, Agudelo-Romero 15, Song 22).",
                       length(unique(gs_lit$gene_set))),
    theme = theme(plot.title    = element_text(family = "Arial", face = "bold",
                                               size = 12, hjust = 0.5,
                                               margin = margin(b = 2)),
                  plot.subtitle = element_text(family = "Arial", size = 9,
                                               hjust = 0.5, colour = "grey25",
                                               margin = margin(b = 6)))
  )

save_grbv(layout, OUT_STEM,
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 200)

cat(sprintf("[S70] %d literature-curated sets across %d categories\n",
            length(unique(gs_lit$gene_set)), nrow(cat_count)))
print(cat_count)
