# figures/scripts_R/S93_hsf_mbf1c_shsp_axis.R
# S93: Custom flow -- HSF -> MBF1c -> sHSP axis (passive-failure model).
# Three-column ggplot layout:
#   Left   : HSF transcription factors from heat_stress_DEGs_all_conditions.csv
#   Middle : MBF1c coactivator (LOC100249249)
#   Right  : sHSP target genes from the same CSV
# Each node carries an 8-cell mini-bar of log2FC across Q2x (pure stage order:
#   2021 Pre / Ver / Post / Harv | 2022 Pre / Ver / Post / Harv).
# Edges: HSFs -> MBF1c, plus HSFs + MBF1c -> all sHSPs.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(ggrepel)
})
source("figures/scripts_R/theme_grbv.R")
source("figures/scripts_R/pathway_overlay_driver.R")

# 1. Q2x log2FC matrix (rownames are stripped Entrez IDs)
M <- read_q2x_matrix()

# 2. Curated gene lists from heat_stress_DEGs CSV
hs <- read_csv("heat_stress_DEGs_all_conditions.csv", show_col_types = FALSE)

# HSFs: 13 unique LOCs in Category == "HSF transcription factors"
hsf_df <- hs %>%
  filter(Category == "HSF transcription factors") %>%
  distinct(Gene, Gene_Name) %>%
  filter(grepl("^LOC\\d+$", Gene))

# Construct short labels from Gene_Name
# e.g. "heat stress transcription factor A-2e" -> "HsfA-2e"
#      "heat shock factor protein HSF24"       -> "HSF24"
#      "heat shock factor-binding protein"     -> "HSF-binding"
make_hsf_label <- function(nm) {
  out <- nm
  out <- gsub("^heat stress transcription factor ?", "Hsf", out, ignore.case = TRUE)
  out <- gsub("^heat shock factor protein ?", "", out, ignore.case = TRUE)
  out <- gsub("^heat shock factor-binding protein.*$", "HSF-binding", out, ignore.case = TRUE)
  out <- substr(out, 1, 18)
  out
}
hsf_df$label  <- vapply(hsf_df$Gene_Name, make_hsf_label, character(1))
# Disambiguate duplicate HSF-binding labels by suffixing the trailing LOC digits.
dup_mask <- duplicated(hsf_df$label) | duplicated(hsf_df$label, fromLast = TRUE)
hsf_df$label[dup_mask] <- paste0(hsf_df$label[dup_mask], " (", substr(hsf_df$Gene[dup_mask], 7, 12), ")")
hsf_df$loc    <- hsf_df$Gene
hsf_df$entrez <- sub("^LOC", "", hsf_df$loc)
hsf_df <- hsf_df[hsf_df$entrez %in% rownames(M), ]
n_hsf <- nrow(hsf_df)
stopifnot("expected >=10 HSFs in matrix" = n_hsf >= 10)

# MBF1c: hardcoded (not in the heat_stress CSV categories)
mbf_df <- tibble(label  = "MBF1c\n(LOC100249249)",
                 loc    = "LOC100249249",
                 entrez = "100249249")
stopifnot("MBF1c not in matrix" = "100249249" %in% rownames(M))

# sHSPs: filter Category == "sHSP (small heat shock proteins)"
shsp_df <- hs %>%
  filter(Category == "sHSP (small heat shock proteins)") %>%
  distinct(Gene, Gene_Name) %>%
  filter(grepl("^LOC\\d+$", Gene))

# Label cleanup: drop trailing organelle qualifiers; collapse "small " prefix to "s"
make_shsp_label <- function(nm) {
  out <- nm
  out <- gsub("\\n", " ", out, fixed = TRUE)
  out <- gsub(",.*$", "", out)            # drop ", chloroplastic-like" etc.
  out <- gsub("^small ", "s", out, ignore.case = TRUE)
  out <- gsub('^"', "", out)              # strip leading quote leftover from CSV oddities
  out <- gsub("kDa heat shock protein", "kDa HSP", out)
  out <- substr(out, 1, 28)
  out
}
shsp_df$label  <- vapply(shsp_df$Gene_Name, make_shsp_label, character(1))
# Disambiguate duplicate labels with trailing LOC digits
dup_mask <- duplicated(shsp_df$label) | duplicated(shsp_df$label, fromLast = TRUE)
shsp_df$label[dup_mask] <- paste0(shsp_df$label[dup_mask], " (", substr(shsp_df$Gene[dup_mask], 7, 12), ")")
shsp_df$loc    <- shsp_df$Gene
shsp_df$entrez <- sub("^LOC", "", shsp_df$loc)
shsp_df <- shsp_df[shsp_df$entrez %in% rownames(M), ]
n_shsp <- nrow(shsp_df)
stopifnot("expected >=20 sHSPs in matrix" = n_shsp >= 20)

# 3. Node tibble with x,y coordinates
nodes <- bind_rows(
  tibble(label = hsf_df$label, loc = hsf_df$loc, entrez = hsf_df$entrez,
         x = 0, y = seq(-n_hsf/2,  n_hsf/2,  length.out = n_hsf),
         group = "HSF"),
  tibble(label = mbf_df$label, loc = mbf_df$loc, entrez = mbf_df$entrez,
         x = 3, y = 0, group = "Coactivator"),
  tibble(label = shsp_df$label, loc = shsp_df$loc, entrez = shsp_df$entrez,
         x = 6, y = seq(-n_shsp/2, n_shsp/2, length.out = n_shsp),
         group = "sHSP")
)

# 4. Edge list: HSFs -> MBF1c; HSFs + MBF1c -> sHSPs
edges <- bind_rows(
  expand.grid(from = hsf_df$label, to = mbf_df$label, stringsAsFactors = FALSE),
  expand.grid(from = c(hsf_df$label, mbf_df$label),
              to   = shsp_df$label, stringsAsFactors = FALSE)
)
edges <- edges %>%
  left_join(nodes %>% select(label, x_from = x, y_from = y), by = c("from" = "label")) %>%
  left_join(nodes %>% select(label, x_to   = x, y_to   = y), by = c("to"   = "label"))

# 5. Per-node 8-cell mini-bar polygons
make_cells <- function(node_row, M, cell_w = 0.7, cell_h = 0.30) {
  ez <- node_row$entrez
  if (!ez %in% rownames(M)) return(NULL)
  fc <- M[ez, ]
  x0 <- node_row$x - cell_w/2
  tibble(
    cell_idx  = 1:8,
    condition = colnames(M),
    fc        = as.numeric(fc),
    xmin      = x0 + (0:7) * (cell_w / 8),
    xmax      = x0 + (1:8) * (cell_w / 8),
    ymin      = node_row$y - cell_h/2,
    ymax      = node_row$y + cell_h/2,
    label     = node_row$label
  )
}
cells <- bind_rows(lapply(seq_len(nrow(nodes)),
                          function(i) make_cells(nodes[i, ], M)))
cells$fc_clip <- pmax(-3, pmin(3, cells$fc))

# 6. Render
p <- ggplot() +
  # Edges first (under nodes), low alpha so dense column-to-column lines don't overpower
  geom_segment(data = edges,
               aes(x = x_from + 0.4, xend = x_to - 0.4,
                   y = y_from,       yend = y_to),
               colour = "grey80", linewidth = 0.12, alpha = 0.4) +
  # 8-cell mini-bars per node
  geom_rect(data = cells,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fc_clip),
            colour = NA) +
  # Node frame
  geom_rect(data = nodes,
            aes(xmin = x - 0.35, xmax = x + 0.35,
                ymin = y - 0.155, ymax = y + 0.155),
            fill = NA, colour = "grey20", linewidth = 0.25) +
  # Labels: HSFs to the left, sHSPs to the right, MBF1c centered above
  geom_text(data = nodes,
            aes(x = x + ifelse(group == "HSF", -0.4,
                        ifelse(group == "sHSP", 0.4, 0)),
                y = y + ifelse(group == "Coactivator", 0.30, 0),
                label = label,
                hjust = ifelse(group == "HSF", 1,
                        ifelse(group == "sHSP", 0, 0.5))),
            size = 2.2, family = "Source Sans Pro") +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                       midpoint = 0, limits = c(-3, 3),
                       oob = scales::squish,
                       name = "log2 FC") +
  coord_cartesian(clip = "off", xlim = c(-2, 8)) +
  labs(title = "S93. HSF -> MBF1c -> sHSP axis",
       subtitle = "Each node: 8-cell log2FC mini-bar (2021 Pre/Ver/Post/Harv | 2022 Pre/Ver/Post/Harv). Edges: HSFs->MBF1c, HSFs+MBF1c->sHSPs.") +
  theme_grbv() +
  theme(axis.text   = element_blank(), axis.ticks = element_blank(),
        axis.title  = element_blank(),
        panel.grid  = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(8, 30, 8, 30))

# 7. Save
ggsave("figures/supp/S93_hsf_mbf1c_shsp_axis.png", p,
       width = 9, height = 11, dpi = 600, units = "in")
ggsave("figures/supp/S93_hsf_mbf1c_shsp_axis.svg", p,
       width = 9, height = 11, units = "in")

cat(sprintf("S93 done: figures/supp/S93_hsf_mbf1c_shsp_axis.png  (%d HSF + %d sHSP nodes)\n",
            n_hsf, n_shsp))
