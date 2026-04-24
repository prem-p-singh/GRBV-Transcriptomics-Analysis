# S69 — Pathway enrichment across infection cohort year (RB+ vs RB-).
#
# Companion to S68 (season x stage enrichment). Answers the year axis.
# For each of 5 cohorts (2016, 2019, 2020, 2021, 2022), we run Fisher-exact
# enrichment of that cohort's DEGs against 13 curated gene-set families:
#
#   Heat-stress categories (9): sHSP, HSP70, HSP90, HSP100/ClpB, HSF,
#                               antioxidant, LEA, GolS, GO:heat
#   TF families (3):            WRKY, NAC, bHLH (HSF already covered above)
#   Infection signature (1):    Tier-1 GRBV 55-gene signature
#
# DEGs per cohort: padj<0.05, |log2FC|>=1 from Q1b_cohort*_deseq2.csv.
# Background universe: all genes tested in that cohort (|log2FC|>0 with padj).
#
# Output: figures/supp/S69_enrichment_cohort.png + .svg

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table); library(stringr)
})

HERE     <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ROOT     <- "/Users/prempratapsingh/Data_Analysis"
FARM     <- file.path(ROOT, "farm_pipeline/results/de_extras")
OUT_DIR  <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S69_enrichment_cohort")

COHORTS <- c("2016","2019","2020","2021","2022")
PADJ    <- 0.05; LFC <- 1.0

# ============================================================================
# LOAD gene-set families
# ============================================================================
heat <- fread(file.path(ROOT, "heat_stress_DEGs_all_conditions.csv")) %>%
  as.data.frame()
# Simplify category names for display
heat$Family <- str_replace_all(
  heat$Category,
  c("small heat shock proteins"   = "sHSP",
    "HSP70 family"                = "HSP70",
    "HSP90 family"                = "HSP90",
    "HSP100/ClpB"                 = "HSP100/ClpB",
    "HSF transcription factors"   = "HSF",
    "Antioxidant enzymes"         = "Antioxidant",
    "LEA/Dehydrin"                = "LEA/dehydrin",
    "Galactinol synthase"         = "GolS",
    "GO:0009408 \\(response to heat\\)" = "GO:heat"))
heat$Family <- gsub("\\(.+\\)|\\s+$|\\s*\\u2022.*$", "", heat$Family) %>% trimws()
heat_sets <- heat %>% select(Gene, Family) %>% distinct() %>%
  filter(Family != "")

# TF families (from existing per-condition TF table)
tf_fam <- fread(file.path(ROOT, "tf_family_gene_details.csv")) %>% as.data.frame()
tf_sets <- tf_fam %>%
  filter(tf_family %in% c("WRKY","NAC","bHLH")) %>%
  select(Gene = gene, Family = tf_family) %>% distinct()

# Tier-1 GRBV signature (padj<0.05, |log2FC|>=1)
sig <- read_csv(file.path(ROOT, "tier1_grbv_signature.csv"),
                show_col_types = FALSE) %>%
  filter(!is.na(padj), padj < PADJ, abs(log2FC) >= LFC)
sig_sets <- data.frame(Gene = sig$LOC_id, Family = "GRBV signature")

all_sets <- bind_rows(heat_sets, tf_sets, sig_sets)

# Preserve a tidy family order for the heatmap row axis
family_order <- c("sHSP","HSP70","HSP90","HSP100/ClpB","HSF",
                  "Antioxidant","LEA/dehydrin","GolS","GO:heat",
                  "WRKY","NAC","bHLH","GRBV signature")
family_order <- family_order[family_order %in% unique(all_sets$Family)]

# ============================================================================
# Run enrichment per cohort
# ============================================================================
enrich_cohort <- function(cohort) {
  f <- file.path(FARM, sprintf("Q1b_cohort%s_vs_neg_deseq2.csv", cohort))
  d <- fread(f) %>% as.data.frame()
  d$is_deg <- !is.na(d$padj) & d$padj < PADJ & abs(d$log2FC) >= LFC
  universe <- d$gene                                # all tested
  deg_set  <- d$gene[d$is_deg]
  n_u <- length(universe); n_d <- length(deg_set)

  # Fisher exact per family: 2x2 (DEG vs not) x (in family vs not),
  # restricted to the cohort's universe.
  rows <- lapply(family_order, function(fam) {
    set_g <- all_sets$Gene[all_sets$Family == fam]
    set_in_u <- intersect(set_g, universe)
    a <- sum(deg_set %in% set_in_u)                 # DEG & in family
    b <- length(deg_set) - a                        # DEG & not in family
    c <- length(setdiff(set_in_u, deg_set))         # not DEG & in family
    dd <- n_u - a - b - c                           # not DEG & not in family
    m  <- matrix(c(a, b, c, dd), nrow = 2, byrow = TRUE)
    if (any(m < 0) || sum(m[, 1]) == 0 || sum(m[1, ]) == 0) {
      return(data.frame(cohort = cohort, family = fam, a = a,
                        set_n = length(set_in_u), deg_n = n_d,
                        fold = NA, pvalue = NA, neglog10q = NA))
    }
    fe <- suppressWarnings(fisher.test(m, alternative = "greater"))
    # Fold enrichment: (a/deg_n) / (set_n/n_u)
    exp_f <- length(set_in_u) / n_u
    obs_f <- a / n_d
    fold  <- ifelse(exp_f > 0, obs_f / exp_f, NA)
    data.frame(cohort = cohort, family = fam, a = a,
               set_n = length(set_in_u), deg_n = n_d,
               fold = fold, pvalue = fe$p.value,
               neglog10q = NA_real_)
  })
  do.call(rbind, rows)
}

res <- bind_rows(lapply(COHORTS, enrich_cohort))

# BH-adjust per family across cohorts, then across whole panel
res$qvalue <- p.adjust(res$pvalue, method = "BH")
res$neglog10q <- -log10(pmin(res$qvalue, 1))
res$neglog10q[!is.finite(res$neglog10q)] <- 0
res$family <- factor(res$family, levels = rev(family_order))
res$cohort <- factor(res$cohort, levels = COHORTS)

# ============================================================================
# Panel A — Heatmap of -log10(q) per cohort x family, cell text = fold
# ============================================================================
res$cell_label <- ifelse(!is.na(res$fold) & res$a > 0,
                         sprintf("%.1fx (%d)", res$fold, res$a), "")

max_nl <- max(res$neglog10q, na.rm = TRUE)

pA <- ggplot(res, aes(x = cohort, y = family, fill = neglog10q)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = cell_label),
            family = "Arial", size = 2.3, colour = "grey15") +
  scale_fill_gradient(low = "#F4EDE2",
                      high = unname(GRBV_COLORS["RB_pos"]),
                      name = "-log\u2081\u2080(q)",
                      limits = c(0, max(1, max_nl)), oob = scales::squish) +
  labs(x = "Infection cohort (year first infected)", y = NULL,
       title = "A  Pathway enrichment in DEGs per cohort",
       subtitle = "Fisher-exact one-sided; cell text = fold enrichment (n overlap). Colour = -log\u2081\u2080(BH q).") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        axis.text.x   = element_text(size = 9),
        axis.text.y   = element_text(size = 8),
        legend.key.width  = unit(8, "pt"),
        legend.key.height = unit(14, "pt"))

# ============================================================================
# Panel B — DEG count trajectory across cohorts (context for A)
# ============================================================================
deg_counts <- bind_rows(lapply(COHORTS, function(c) {
  f <- file.path(FARM, sprintf("Q1b_cohort%s_vs_neg_deseq2.csv", c))
  d <- fread(f) %>% as.data.frame()
  d$is_deg <- !is.na(d$padj) & d$padj < PADJ & abs(d$log2FC) >= LFC
  data.frame(cohort = c,
             up   = sum(d$is_deg & d$log2FC > 0),
             down = sum(d$is_deg & d$log2FC < 0))
})) %>%
  pivot_longer(c(up, down), names_to = "dir", values_to = "n") %>%
  mutate(cohort = factor(cohort, levels = COHORTS),
         dir    = factor(dir, levels = c("up","down")))

totals <- deg_counts %>% group_by(cohort) %>%
  summarise(total = sum(n), .groups = "drop")

pB <- ggplot(deg_counts, aes(x = cohort, y = n, fill = dir)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.65,
           colour = "white", linewidth = 0.25) +
  geom_text(data = totals,
            aes(x = cohort, y = total, label = scales::comma(total)),
            inherit.aes = FALSE, vjust = -0.4,
            family = "Arial", size = 2.6) +
  scale_fill_manual(values = c(up = unname(GRBV_COLORS["up"]),
                               down = unname(GRBV_COLORS["down"])),
                    labels = c(up = "Up", down = "Down"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = scales::label_comma()) +
  labs(x = "Cohort", y = "DEGs (padj<0.05, |log\u2082FC|\u22651)",
       title = "B  DEGs per cohort",
       subtitle = "recent cohorts = more DEGs") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        legend.position = "top",
        legend.text   = element_text(size = 8),
        legend.key.height = unit(8, "pt"))

# ============================================================================
# Compose
# ============================================================================
layout <- (pA | pB) + plot_layout(widths = c(1.6, 0.9)) +
  plot_annotation(
    title = "Fig. S69 \u2014 Pathway enrichment shifts with years since infection",
    subtitle = "Fisher-exact enrichment of per-cohort DEGs against 13 curated gene-set families.",
    theme = theme(plot.title    = element_text(family = "Arial", face = "bold",
                                               size = 12, hjust = 0.5,
                                               margin = margin(b = 2)),
                  plot.subtitle = element_text(family = "Arial", size = 9,
                                               hjust = 0.5, colour = "grey25",
                                               margin = margin(b = 6)))
  )

save_grbv(layout, OUT_STEM,
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 150)

cat(sprintf("[S69] %d families x %d cohorts; %d significant (q<0.05) cells\n",
            length(family_order), length(COHORTS),
            sum(res$qvalue < 0.05, na.rm = TRUE)))
