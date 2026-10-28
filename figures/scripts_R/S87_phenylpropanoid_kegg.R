# figures/scripts_R/S87_phenylpropanoid_kegg.R
# S87: Phenylpropanoid biosynthesis (KEGG vvi00940) with 8-cell split-box overlay.
# Direct paper-relevance: anthocyanin / colour-loss phenotype; GSEA 4/8 cond, NES -1.14.

source("figures/scripts_R/pathway_overlay_driver.R")

M <- read_q2x_matrix()
res <- render_kegg_pathway(pathway_id = "vvi00940",
                           supp_id    = 87,
                           slug       = "phenylpropanoid_kegg",
                           M          = M)
cat(sprintf("S87 done: %s (coverage %.1f%%)\n", res$png_path, res$coverage_pct))
