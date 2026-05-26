# figures/scripts_R/S89_plant_hormone_kegg.R
# S89: Plant hormone signal transduction (KEGG vvi04075) with 8-cell split-box overlay.
# Direct paper-relevance: GA suppression (4/4 sig DOWN) and ethylene cross-paper match vs BU 2017.
# GSEA 2/8 cond, NES +1.48 (narrative pathway).

source("figures/scripts_R/pathway_overlay_driver.R")

M <- read_q2x_matrix()
res <- render_kegg_pathway(pathway_id = "vvi04075",
                           supp_id    = 89,
                           slug       = "plant_hormone_kegg",
                           M          = M)
cat(sprintf("S89 done: %s (coverage %.1f%%)\n", res$png_path, res$coverage_pct))
