# figures/scripts_R/S88_starch_sucrose_kegg.R
# S88: Starch and sucrose metabolism (KEGG vvi00500) with 8-cell split-box overlay.
# Direct paper-relevance: harvest Brix suppression and SPS1 mediator reframe.
# GSEA 5/8 cond, NES +0.89 (UP under GRBV).

source("figures/scripts_R/pathway_overlay_driver.R")

M <- read_q2x_matrix()
res <- render_kegg_pathway(pathway_id = "vvi00500",
                           supp_id    = 88,
                           slug       = "starch_sucrose_kegg",
                           M          = M)
cat(sprintf("S88 done: %s (coverage %.1f%%)\n", res$png_path, res$coverage_pct))
