# figures/scripts_R/S92_proteasome_kegg.R
# S92: Proteasome (KEGG vvi03050) with 8-cell split-box overlay.
# Direct paper-relevance: strongest |NES| in the GSEA panel; ubiquitin-proteasome system collapse.
# GSEA 4/8 cond, NES -2.43 (strongly DOWN under GRBV).

source("figures/scripts_R/pathway_overlay_driver.R")

M <- read_q2x_matrix()
res <- render_kegg_pathway(pathway_id = "vvi03050",
                           supp_id    = 92,
                           slug       = "proteasome_kegg",
                           M          = M)
cat(sprintf("S92 done: %s (coverage %.1f%%)\n", res$png_path, res$coverage_pct))
