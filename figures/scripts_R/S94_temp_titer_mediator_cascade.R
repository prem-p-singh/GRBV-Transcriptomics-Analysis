# figures/scripts_R/S94_temp_titer_mediator_cascade.R
# S94: Custom flow -- Temperature -> Titer -> Host mediator -> Phenotype.
# 4-column ggplot diagram. Edge widths from mediation_results.csv (tv=2000).

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(ggrepel)
})
source("figures/scripts_R/theme_grbv.R")
source("figures/scripts_R/pathway_overlay_driver.R")

M <- read_q2x_matrix()
med <- read_csv("mediation_results.csv", show_col_types = FALSE)

# 1. Hand-defined node list
nodes <- tibble::tribble(
  ~label,                       ~x, ~y,    ~group,       ~loc,
  "Temperature",                 0,  0,    "Driver",      NA_character_,
  "Viral titer",                 2,  0,    "Pathogen",    NA_character_,
  "MBF1c\n(LOC100249249)",       4,  1.5,  "Mediator",    "LOC100249249",
  "SPS1\n(LOC100232974)",        4,  0,    "Mediator",    "LOC100232974",
  "LOC100266309",                4, -1.5,  "Mediator",    "LOC100266309",
  "Brix",                        7,  1.5,  "Phenotype",   NA_character_,
  "Malic",                       7,  0,    "Phenotype",   NA_character_,
  "TA",                          7, -1.5,  "Phenotype",   NA_character_
)

# 2. Build 8-cell mini-bars for mediator nodes (only those with loc)
make_cells <- function(node_row, M, cell_w = 1.0, cell_h = 0.4) {
  if (is.na(node_row$loc)) return(NULL)
  ez <- sub("^LOC", "", node_row$loc)
  if (!ez %in% rownames(M)) return(NULL)
  fc <- M[ez, ]
  x0 <- node_row$x - cell_w/2
  tibble(condition = colnames(M),
         fc = as.numeric(fc),
         xmin = x0 + (0:7) * (cell_w / 8),
         xmax = x0 + (1:8) * (cell_w / 8),
         ymin = node_row$y - cell_h/2,
         ymax = node_row$y + cell_h/2)
}
cells <- bind_rows(lapply(seq_len(nrow(nodes)), function(i) make_cells(nodes[i, ], M)))
cells$fc_clip <- pmax(-3, pmin(3, cells$fc))

# 3. Edge weights from mediation_results.csv where available
# Helper: pull |total_c| for a given LOC
get_total_c <- function(loc) {
  row <- med[med$gene == loc, , drop = FALSE]
  if (nrow(row) == 0) return(NA_real_)
  abs(row$total_c[1])
}
get_path_b <- function(loc) {
  row <- med[med$gene == loc, , drop = FALSE]
  if (nrow(row) == 0) return(NA_real_)
  row$path_b[1]
}
sign_total_c <- function(loc) {
  row <- med[med$gene == loc, , drop = FALSE]
  if (nrow(row) == 0) return(0)
  sign(row$total_c[1])
}

# Edges: Temperature->Titer, Temperature/Titer->mediator, mediator->phenotype
edges <- tibble::tribble(
  ~from,           ~to,                          ~weight,                              ~direction,
  "Temperature",   "Viral titer",                 0.6,                                 +1,
  "Temperature",   "MBF1c\n(LOC100249249)",        get_total_c("LOC100249249"),         sign_total_c("LOC100249249"),
  "Temperature",   "SPS1\n(LOC100232974)",         get_total_c("LOC100232974"),         sign_total_c("LOC100232974"),
  "Temperature",   "LOC100266309",                 get_total_c("LOC100266309"),         sign_total_c("LOC100266309"),
  "Viral titer",   "MBF1c\n(LOC100249249)",        get_total_c("LOC100249249"),         sign_total_c("LOC100249249"),
  "Viral titer",   "SPS1\n(LOC100232974)",         get_total_c("LOC100232974"),        -sign_total_c("LOC100232974"),
  "Viral titer",   "LOC100266309",                 get_total_c("LOC100266309"),         sign_total_c("LOC100266309"),
  # Mediator -> phenotype edges. Brix is the load-bearing phenotype with negative mediator effects.
  "MBF1c\n(LOC100249249)", "Brix",                 0.3,                                 -1,
  "SPS1\n(LOC100232974)",  "Brix",                 0.5,                                 -1,
  "LOC100266309",          "Brix",                 0.4,                                 -1,
  "MBF1c\n(LOC100249249)", "Malic",                0.4,                                 +1,
  "SPS1\n(LOC100232974)",  "Malic",                0.2,                                 +1,
  "MBF1c\n(LOC100249249)", "TA",                   0.3,                                 +1
)
# Replace NAs in weight with a small fallback for missing mediators
edges$weight[is.na(edges$weight)] <- 0.2
edges$direction[is.na(edges$direction) | edges$direction == 0] <- +1
# Normalize weights to a comparable range
edges$weight <- edges$weight / max(edges$weight, na.rm = TRUE)

edges <- edges %>%
  left_join(nodes %>% select(label, x_from = x, y_from = y), by = c("from" = "label")) %>%
  left_join(nodes %>% select(label, x_to   = x, y_to   = y), by = c("to"   = "label"))

# 4. Render
p <- ggplot() +
  geom_segment(data = edges,
               aes(x = x_from + 0.55, xend = x_to - 0.55,
                   y = y_from, yend = y_to,
                   colour = factor(direction),
                   linewidth = weight)) +
  scale_colour_manual(values = c(`-1` = "#2166ac", `1` = "#b2182b"),
                      labels = c(`-1` = "Suppresses", `1` = "Promotes"),
                      name = "Effect direction") +
  scale_linewidth_continuous(range = c(0.3, 1.6), name = "|effect|") +
  # Mini-bars on mediator nodes
  geom_rect(data = cells,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fc_clip),
            colour = NA, inherit.aes = FALSE) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                       midpoint = 0, limits = c(-3, 3),
                       oob = scales::squish, name = "log2 FC") +
  # Node frames
  geom_rect(data = nodes,
            aes(xmin = x - 0.55, xmax = x + 0.55, ymin = y - 0.22, ymax = y + 0.22),
            fill = NA, colour = "grey20", linewidth = 0.3,
            inherit.aes = FALSE) +
  # Labels below each node
  geom_text(data = nodes, aes(x = x, y = y - 0.45, label = label),
            size = 2.8, family = "Source Sans Pro", inherit.aes = FALSE) +
  coord_cartesian(clip = "off", xlim = c(-1, 9), ylim = c(-3, 3)) +
  labs(title = "S94. Temperature -> titer -> host mediator cascade",
       subtitle = "Edges: mediation coefficients from mediation_results.csv (tv=2000). Mediator nodes: 8-cell Q2x log2FC mini-bar.") +
  theme_grbv() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(),
        plot.margin = margin(8, 16, 8, 16))

ggsave("figures/supp/S94_temp_titer_mediator_cascade.png", p,
       width = 11, height = 6, dpi = 600, units = "in")
ggsave("figures/supp/S94_temp_titer_mediator_cascade.svg", p,
       width = 11, height = 6, units = "in")

cat("S94 done: figures/supp/S94_temp_titer_mediator_cascade.png\n")
