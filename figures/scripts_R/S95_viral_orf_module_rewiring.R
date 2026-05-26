# figures/scripts_R/S95_viral_orf_module_rewiring.R
# S95: Custom bipartite -- 6 viral ORFs -> 13 WGCNA modules.
# Edges weighted by per-ORF x per-module Spearman (rule 13). Filter: padj < 0.05 OR |rho| >= 0.3.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(ggrepel)
})
source("figures/scripts_R/theme_grbv.R")

# 1. Load ORF x module Spearman correlations from rule 13
orf_mod <- read_csv("farm_pipeline/results/13_viral_orf_ext/orf_module_correlation.csv",
                    show_col_types = FALSE)

# Strip ME prefix from module column for display
orf_mod$module_short <- sub("^ME", "", orf_mod$module)

# 2. ORF nodes (left column)
orfs <- sort(unique(orf_mod$ORF))
orf_nodes <- tibble(label = orfs,
                    short = sub("^N761_gp", "V", orfs),
                    x = 0,
                    y = seq(-length(orfs)/2, length(orfs)/2, length.out = length(orfs)))

# 3. Module nodes (right column)
mods <- sort(unique(orf_mod$module_short))
mod_nodes <- tibble(label = mods,
                    x = 5,
                    y = seq(-length(mods)/2, length(mods)/2, length.out = length(mods)))

# 4. Edges: keep significant or strong correlations
edges <- orf_mod %>%
  filter(padj < 0.05 | abs(rho) >= 0.3) %>%
  left_join(orf_nodes %>% select(orf = label, x_from = x, y_from = y), by = c("ORF" = "orf")) %>%
  left_join(mod_nodes %>% select(module_short = label, x_to = x, y_to = y), by = "module_short") %>%
  mutate(direction = sign(rho), weight = abs(rho))

# 5. Map WGCNA module names to actual colours for the right-column dots
# Most module names ARE valid R colours; "greenyellow" is unusual but valid.
mod_colours <- setNames(mods, mods)  # use the module name as the colour name directly
# Fall back to grey for any name R doesn't recognize
valid_colours <- vapply(mods, function(m) tryCatch({ col2rgb(m); TRUE }, error = function(e) FALSE),
                        logical(1))
mod_colours[!valid_colours] <- "grey50"
names(mod_colours) <- mods

# 6. Render
p <- ggplot() +
  geom_segment(data = edges,
               aes(x = x_from + 0.2, xend = x_to - 0.2,
                   y = y_from, yend = y_to,
                   colour = factor(direction),
                   linewidth = weight),
               alpha = 0.75) +
  scale_colour_manual(values = c(`-1` = "#2166ac", `1` = "#b2182b"),
                      labels = c(`-1` = "Negative", `1` = "Positive"),
                      name = "Spearman rho") +
  scale_linewidth_continuous(range = c(0.25, 1.4), name = "|rho|") +
  # ORF nodes: navy dots on left
  geom_point(data = orf_nodes, aes(x = x, y = y), size = 5.5, colour = "#1f2937") +
  geom_text(data = orf_nodes, aes(x = x - 0.25, y = y, label = short),
            hjust = 1, size = 3.2, family = "Source Sans Pro") +
  # Module nodes: coloured by module name on right
  geom_point(data = mod_nodes, aes(x = x, y = y), size = 5,
             fill = mod_colours[mod_nodes$label], shape = 21, stroke = 0.3) +
  geom_text(data = mod_nodes, aes(x = x + 0.25, y = y, label = label),
            hjust = 0, size = 3.2, family = "Source Sans Pro") +
  coord_cartesian(clip = "off", xlim = c(-1.5, 7.5)) +
  labs(title = "S95. Viral ORFs -> WGCNA module rewiring",
       subtitle = "Per-ORF x per-module Spearman rho (rule 13). Edges shown where FDR < 0.05 or |rho| >= 0.3.") +
  theme_grbv() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(),
        plot.margin = margin(8, 16, 8, 16))

ggsave("figures/supp/S95_viral_orf_module_rewiring.png", p,
       width = 9, height = 8, dpi = 600, units = "in")
ggsave("figures/supp/S95_viral_orf_module_rewiring.svg", p,
       width = 9, height = 8, units = "in")

cat(sprintf("S95 done: figures/supp/S95_viral_orf_module_rewiring.png  (%d ORFs x %d modules, %d edges drawn)\n",
            length(orfs), length(mods), nrow(edges)))
