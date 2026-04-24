# Supplementary Tier 1 literature meta-analysis panels.
#
# Reads the CSVs published by rule 12 (farm_pipeline/scripts/12_comparative_metaanalysis.py)
# at the repo root (tier1_*.csv) and produces two supplementary figures:
#
#   S64 — Tier 1 DEG-overlap summary (3 panels)
#         Panel A: signature size per dataset (log-x horizontal bars)
#         Panel B: GRBV-signature intersections with each comparator
#         Panel C: Fisher-exact enrichment of GRBV vs each comparator
#
#   S65 — Tier 1 log2-fold-change concordance heatmap
#         Square Spearman-rho heatmap across all pairs of Tier 1 datasets.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(ggrepel)
  library(stringr)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

DATA_ROOT <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR   <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- load ----------------------------------------------------------------
overlap  <- read_csv(file.path(DATA_ROOT, "tier1_overlap_matrix.csv"),
                     show_col_types = FALSE)
coverage <- read_csv(file.path(DATA_ROOT, "tier1_coverage_report.csv"),
                     show_col_types = FALSE)
fc_mat   <- read_csv(file.path(DATA_ROOT, "tier1_fc_correlation_matrix.csv"),
                     show_col_types = FALSE) %>%
              rename(dataset = 1) %>% as.data.frame()
spec_cls <- read_csv(file.path(DATA_ROOT, "tier1_specificity_classes.csv"),
                     show_col_types = FALSE)

grbv_sig <- read_csv(file.path(DATA_ROOT, "tier1_grbv_signature.csv"),
                     show_col_types = FALSE)
PADJ <- 0.05; LFC <- 1.0
grbv_set <- grbv_sig %>% filter(padj < PADJ, abs(log2FC) >= LFC) %>% pull(LOC_id)

# ---- short display labels -----------------------------------------------
short_label <- function(x) {
  if (x == "GRBV") return("GRBV (this study)")
  m <- str_match(x, "^([A-Za-z]+)(\\d{4})_([^_]+)_([^_]+)_([^_]+)(?:_(.+))?$")
  if (is.na(m[1,1])) return(x)
  auth <- m[1,2]; yr <- m[1,3]; stress <- m[1,4]
  cult <- m[1,5] %>%
    str_replace("CabernetSauvignon", "CS") %>%
    str_replace("CabernetFranc", "CF") %>%
    str_replace("Zinfandel", "Zin") %>%
    str_replace("Trincadeira", "Trin") %>%
    str_replace("Vvinifera", "Vv")
  tissue <- m[1,6]; stage <- m[1,7]
  paste0(auth, " ", yr, "\n", stress, ", ", cult, " ", tissue,
         if (!is.na(stage)) paste0(" ", stage) else "")
}

# =========================================================================
# S64 — overlap summary (3 panels)
# =========================================================================

# Panel A: signature size per dataset
# Reconstruct sig-sizes per dataset from the overlap matrix (sig_A column).
sig_sizes <- overlap %>%
  select(dataset = dataset_A, n = sig_A) %>%
  bind_rows(overlap %>% select(dataset = dataset_B, n = sig_B)) %>%
  distinct() %>%
  group_by(dataset) %>% summarise(n = max(n, na.rm = TRUE), .groups = "drop")

# Add datasets not appearing in overlap_matrix (zero-coverage papers) via
# the coverage report so we show the full 6-paper panel.
all_datasets <- unique(c(sig_sizes$dataset, coverage$dataset, "GRBV"))
sig_sizes <- tibble(dataset = all_datasets) %>%
  left_join(sig_sizes, by = "dataset") %>%
  mutate(n = replace_na(n, 0))

sig_sizes$label <- vapply(sig_sizes$dataset, short_label, character(1))
sig_sizes <- sig_sizes %>% arrange(n)
sig_sizes$label <- factor(sig_sizes$label, levels = sig_sizes$label)

A <- ggplot(sig_sizes, aes(x = pmax(n, 0.5), y = label)) +
  geom_col(fill = "#6C7A89", width = 0.7) +
  geom_text(aes(label = n),
            hjust = -0.15, size = 2.6, colour = "grey25") +
  scale_x_log10(expand = expansion(mult = c(0, 0.35))) +
  labs(title = "Signature size",
       x = "DEGs at padj < 0.05, |log2FC| \u2265 1",
       y = NULL) +
  theme_grbv(base_size = 9) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 8))

# Panel B: GRBV-signature intersections with each comparator.
# Build from specificity_classes table (which records sig columns per dataset).
sig_cols <- grep("_sig$", colnames(spec_cls), value = TRUE)
others <- setdiff(sub("_sig$", "", sig_cols), "GRBV")

grbv_mat <- spec_cls %>% filter(LOC_id %in% grbv_set)
# recompute sig bool per row per comparator
comparator_sig <- lapply(others, function(o) {
  col <- paste0(o, "_sig")
  ifelse(is.na(grbv_mat[[col]]), FALSE, grbv_mat[[col]] == 1)
})
names(comparator_sig) <- others

grbv_only <- sum(!Reduce(`|`, comparator_sig))
pairwise  <- sapply(others, function(o) {
  # exclusive: this comparator only, none of the others
  in_o   <- comparator_sig[[o]]
  others_bool <- comparator_sig[setdiff(others, o)]
  if (length(others_bool)) in_others <- Reduce(`|`, others_bool) else in_others <- rep(FALSE, length(in_o))
  sum(in_o & !in_others)
})
multi_way <- sum(Reduce(`+`, comparator_sig) >= 2)

rowsB <- tibble(
  label = c("GRBV-only",
            vapply(others, function(o) paste0("GRBV \u2229 ", strsplit(short_label(o), "\n")[[1]][1]), character(1)),
            "GRBV \u2229 \u22652 others"),
  n     = c(grbv_only, unname(pairwise), multi_way),
  kind  = c("only", rep("pair", length(others)), "multi")
) %>% filter(n > 0) %>% arrange(n)
rowsB$label <- factor(rowsB$label, levels = rowsB$label)

kind_pal <- c(only = unname(GRBV_COLORS["med_titer"]),
              pair = "#5DC863",
              multi = unname(GRBV_COLORS["RB_pos"]))

B <- ggplot(rowsB, aes(x = n, y = label, fill = kind)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), hjust = -0.2, size = 2.7, colour = "grey25") +
  scale_fill_manual(values = kind_pal, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = sprintf("GRBV intersections (n = %d)", length(grbv_set)),
       x = "Number of GRBV DEGs", y = NULL) +
  theme_grbv(base_size = 9) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 8))

# Panel C: Fisher enrichment of GRBV vs each comparator
grbv_rows <- overlap %>%
  filter(dataset_A == "GRBV" | dataset_B == "GRBV") %>%
  mutate(other = if_else(dataset_A == "GRBV", dataset_B, dataset_A),
         log2_or = log2(pmax(odds_ratio, 1e-3)),
         q_label = sprintf("q = %.2g", fisher_q)) %>%
  arrange(odds_ratio)
grbv_rows$other_lab <- vapply(grbv_rows$other, short_label, character(1))
grbv_rows$other_lab <- factor(grbv_rows$other_lab, levels = grbv_rows$other_lab)

# Panel C uses a wider x-range to keep the q-value labels inside the plot
# area (avoids clipping by the patchwork/panel boundary).
xmax_C <- max(grbv_rows$log2_or, na.rm = TRUE)
C <- ggplot(grbv_rows, aes(x = log2_or, y = other_lab)) +
  geom_col(fill = unname(GRBV_COLORS["med_titer"]), width = 0.65) +
  geom_text(aes(label = q_label, x = log2_or + xmax_C * 0.05),
            hjust = 0, size = 2.6, colour = "grey25") +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
  scale_x_continuous(limits = c(0, xmax_C * 1.55), expand = c(0, 0)) +
  labs(title = "Fisher enrichment",
       x = "log2(odds ratio)", y = NULL) +
  theme_grbv(base_size = 9) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 8),
        plot.margin = margin(4, 10, 4, 6))

S64 <- (A + B + C) +
  plot_layout(ncol = 3, widths = c(1.0, 1.0, 1.0)) +
  plot_annotation(
    tag_levels = "a",
    title = "Fig. S64 \u2014 Tier 1 DEG overlap: signature size, intersections, enrichment",
    subtitle = "6-paper comparative panel vs GRBV (this study); n = 55 GRBV DEGs.",
    theme = theme_grbv(base_size = 10) +
      theme(plot.title = element_text(face = "bold", size = 11),
            plot.subtitle = element_text(size = 9, colour = "grey25"),
            plot.tag = element_text(face = "bold", size = 11))
  )

save_grbv(S64, file.path(OUT_DIR, "S64_tier1_overlap_summary"),
          width_mm = 280, height_mm = 115)

# =========================================================================
# S65 — log2FC concordance heatmap
# =========================================================================
mat <- fc_mat
rownames(mat) <- mat$dataset
mat$dataset <- NULL
# Coerce to numeric matrix
m <- as.matrix(mat)
storage.mode(m) <- "numeric"

# Short labels on row/col
lbl <- vapply(rownames(m), short_label, character(1))
rownames(m) <- lbl
colnames(m) <- lbl

# Long format for ggplot
d <- as.data.frame(as.table(m))
names(d) <- c("rowv", "colv", "rho")

# Keep axis order consistent with input
ord <- lbl
d$rowv <- factor(d$rowv, levels = ord)
d$colv <- factor(d$colv, levels = ord)

# v2 fix: NA cells now show "no overlap" text so readers can tell absence
# from uninformative silence; add a coverage side-strip showing ID-space
# mapping quality per dataset (why some cells are NA).
heat <- ggplot(d, aes(x = colv, y = rowv, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(is.na(rho), "no\noverlap",
                               sprintf("%.2f", rho)),
                colour = is.na(rho)),
            size = 2.3, lineheight = 0.85, family = "Arial") +
  scale_fill_gradient2(low = unname(GRBV_COLORS["down"]),
                       mid = "#F7F7F7",
                       high = unname(GRBV_COLORS["up"]),
                       midpoint = 0, limits = c(-1, 1),
                       na.value = "#E6E6E6",
                       name = "Spearman \u03c1\n(log2FC)") +
  scale_colour_manual(values = c(`TRUE` = "grey45", `FALSE` = "grey10"),
                      guide = "none") +
  coord_fixed() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_grbv(base_size = 9) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(4, 4, 4, 4))

# Coverage side strip: one horizontal bar per dataset (matched row order)
# using tier1_coverage_report.csv coverage_pct.  GRBV itself has coverage 100 %.
cov_df <- coverage %>% select(dataset, coverage_pct) %>%
  bind_rows(tibble(dataset = "GRBV", coverage_pct = 100)) %>%
  mutate(label = vapply(dataset, short_label, character(1))) %>%
  filter(label %in% lbl) %>%
  mutate(label = factor(label, levels = lbl))

cov_bar <- ggplot(cov_df, aes(x = coverage_pct, y = label)) +
  geom_col(fill = "#7B7B7B", width = 0.65, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.0f%%", coverage_pct)),
            hjust = -0.25, size = 2.4, family = "Arial", colour = "grey25") +
  scale_x_continuous(limits = c(0, 160), breaks = c(0, 50, 100),
                     expand = c(0, 0)) +
  labs(title = "ID-space coverage",
       x = "% of source DEGs\nresolved to LOC", y = NULL) +
  theme_grbv(base_size = 9) +
  theme(plot.title  = element_text(face = "bold", size = 9),
        axis.text.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.x = element_text(size = 7),
        axis.title.x= element_text(size = 7.5, lineheight = 1),
        panel.grid  = element_blank(),
        plot.margin = margin(4, 6, 4, 2))

S65 <- heat + cov_bar +
  plot_layout(widths = c(1.0, 0.28)) +
  plot_annotation(
    title    = "Fig. S65 \u2014 Tier 1 log2FC concordance across Vitis stress datasets",
    subtitle = "Spearman \u03c1 on shared genes; 'no overlap' = no ID-space intersection. Side strip shows LOC-mapping coverage.",
    theme = theme_grbv(base_size = 10) +
      theme(plot.title    = element_text(face = "bold", size = 11),
            plot.subtitle = element_text(size = 9, colour = "grey25"))
  )

save_grbv(S65, file.path(OUT_DIR, "S65_tier1_fc_concordance"),
          width_mm = 215, height_mm = 170)

cat("supp_tier1.R: wrote S64 and S65 to", OUT_DIR, "\n")
