# figures/scripts_R/S90_ribosome_kegg.R
# S90: Ribosome (KEGG vvi03010) with 8-cell split-box overlay.
# Direct paper-relevance: most robust DOWN signal across the GSEA panel; broad translation suppression.
# GSEA 6/8 cond, NES -1.33 (strongly DOWN under GRBV).

source("figures/scripts_R/pathway_overlay_driver.R")

M <- read_q2x_matrix()
res <- render_kegg_pathway(pathway_id = "vvi03010",
                           supp_id    = 90,
                           slug       = "ribosome_kegg",
                           M          = M)
cat(sprintf("S90 done: %s (coverage %.1f%%)\n", res$png_path, res$coverage_pct))
