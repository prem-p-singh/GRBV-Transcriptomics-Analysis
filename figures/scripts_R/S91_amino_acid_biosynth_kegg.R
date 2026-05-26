# figures/scripts_R/S91_amino_acid_biosynth_kegg.R
# S91: Biosynthesis of amino acids (KEGG vvi01230) with 8-cell split-box overlay.
# Direct paper-relevance: nitrogen metabolism depletion mirrored at harvest (low NOPA/NH3 in late cohorts).
# GSEA 5/8 cond, NES -1.95 (DOWN under GRBV).

source("figures/scripts_R/pathway_overlay_driver.R")

M <- read_q2x_matrix()
res <- render_kegg_pathway(pathway_id = "vvi01230",
                           supp_id    = 91,
                           slug       = "amino_acid_biosynth_kegg",
                           M          = M)
cat(sprintf("S91 done: %s (coverage %.1f%%)\n", res$png_path, res$coverage_pct))
