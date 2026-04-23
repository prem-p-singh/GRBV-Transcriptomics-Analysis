# Supplementary batch: stratified + direction-resolved pathway enrichment
# per research question (S8 - S14) using DESeq2 fix01 / fix01c contrasts.
#
# 2026-04-22 rewrite:
#   * Every enrichment is now run separately on UP- and DOWN-regulated DEGs
#     so the biological interpretation can be direction-aware.
#   * S8 (Q1): heatmap over 5 infection-onset cohorts (Q1b) x direction,
#     recovering the age-trajectory lost in the previous single-bar version.
#   * S9 (Q2): heatmap over 4 maturity stages x direction.
#   * S10-S13: each is a 2x2 grid showing GO:BP UP / DOWN and KEGG UP / DOWN.
#   * S14: master landscape across Q1-Q6 x direction (12 columns).
#   * Dropped the make.unique(..., sep = " #") label anti-pattern.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table)
  library(clusterProfiler); library(org.Vvinifera.eg.db); library(AnnotationDbi)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/.claude/worktrees/gallant-sinoussi-04711f/farm_pipeline/scripts/figures"
source(file.path(HERE, "theme_grbv.R"))

RESULTS <- "/Users/prempratapsingh/Data_Analysis/farm_pipeline/results"
OUT_DIR <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PADJ  <- 0.05
LFC   <- 1.0
N_TOP <- 12        # terms per direction shown in bar panels
N_MASTER <- 25     # terms shown in master / stratified heatmaps

# -----------------------------------------------------------------
# Load every DESeq2 table
# -----------------------------------------------------------------
all_de <- rbindlist(lapply(c(
  list.files(file.path(RESULTS, "fix01_de_redo"),
             pattern = "_deseq2\\.csv$", full.names = TRUE),
  list.files(file.path(RESULTS, "fix01c_de_extras"),
             pattern = "_deseq2\\.csv$", full.names = TRUE)
), function(f) { d <- fread(f); d$file <- basename(f); d }),
fill = TRUE) %>% as.data.frame()

universe <- unique(all_de$gene)

# org.Vvinifera.eg.db uses Entrez IDs without the "LOC" prefix.
strip_loc <- function(ids) sub("^LOC", "", ids)
universe_entrez <- strip_loc(universe)
entrez_to_loc   <- setNames(universe, universe_entrez)

# -----------------------------------------------------------------
# Build TERM2GENE + TERM2NAME for GO BP
# -----------------------------------------------------------------
message("[supp_pathway] fetching GO BP mappings ...")
go_map <- AnnotationDbi::select(org.Vvinifera.eg.db,
                                keys = universe_entrez, keytype = "GID",
                                columns = c("GO","ONTOLOGY"))
go_map <- go_map[!is.na(go_map$GO) & go_map$ONTOLOGY == "BP", ]
go_map$LOC <- entrez_to_loc[go_map$GID]
term2gene_bp <- go_map[!is.na(go_map$LOC), c("GO","LOC")]
colnames(term2gene_bp) <- c("term","gene")

suppressWarnings(library(GO.db))
go_terms <- AnnotationDbi::select(GO.db,
                                  keys = unique(term2gene_bp$term),
                                  columns = c("GOID","TERM","ONTOLOGY"),
                                  keytype = "GOID")
term2name_bp <- go_terms[, c("GOID","TERM")]
colnames(term2name_bp) <- c("term","name")
message(sprintf("[supp_pathway] %d GO:BP term->gene pairs, %d terms",
                nrow(term2gene_bp), length(unique(term2gene_bp$term))))

# -----------------------------------------------------------------
# KEGG
# -----------------------------------------------------------------
kegg_map <- AnnotationDbi::select(org.Vvinifera.eg.db,
                                  keys = universe_entrez, keytype = "GID",
                                  columns = "KEGG")
kegg_map <- kegg_map[!is.na(kegg_map$KEGG), ]
kegg_map$LOC <- entrez_to_loc[kegg_map$GID]
term2gene_kg <- kegg_map[!is.na(kegg_map$LOC), c("KEGG","LOC")]
colnames(term2gene_kg) <- c("term","gene")
term2name_kg <- data.frame(term = unique(term2gene_kg$term),
                            name = paste0("KEGG:", unique(term2gene_kg$term)))
message(sprintf("[supp_pathway] %d KEGG term->gene pairs", nrow(term2gene_kg)))

# -----------------------------------------------------------------
# Enrichment primitives
# -----------------------------------------------------------------
enrich_one <- function(sig_genes, t2g, t2n) {
  if (length(sig_genes) == 0) return(NULL)
  suppressWarnings({
    e <- enricher(
      gene          = sig_genes,
      universe      = universe,
      TERM2GENE     = t2g,
      TERM2NAME     = t2n,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.20,
      minGSSize     = 5,
      maxGSSize     = 500
    )
  })
  if (is.null(e) || nrow(as.data.frame(e)) == 0) return(NULL)
  as.data.frame(e) %>%
    dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count) %>%
    # Some enricher rows have Description as NA or empty when TERM2NAME
    # lookup misses (typically obsolete GO IDs). Drop them at the source.
    filter(!is.na(Description), Description != "", Description != "NA",
           nchar(trimws(Description)) > 0)
}

# Split a DE table into UP / DOWN gene lists
split_dir <- function(d) {
  ok <- !is.na(d$padj) & d$padj < PADJ & abs(d$log2FC) >= LFC
  list(UP   = d$gene[ok & d$log2FC >=  LFC],
       DOWN = d$gene[ok & d$log2FC <= -LFC])
}

# Given a single DE table (one comparison), return direction-resolved
# enrichment for both ontologies as a tall data frame.
enrich_by_direction <- function(d) {
  sets <- split_dir(d)
  rows <- list()
  for (dir in c("UP","DOWN")) {
    for (onto in c("GO:BP","KEGG")) {
      t2g <- if (onto == "GO:BP") term2gene_bp else term2gene_kg
      t2n <- if (onto == "GO:BP") term2name_bp else term2name_kg
      r <- enrich_one(sets[[dir]], t2g, t2n)
      if (!is.null(r) && nrow(r))
        rows[[length(rows)+1]] <- r %>%
          mutate(direction = dir, ontology = onto)
    }
  }
  if (length(rows) == 0) return(NULL)
  bind_rows(rows)
}

# ------------ pick DE rows for a named contrast --------------
de_of <- function(comp) all_de %>% filter(comparison == comp)

# =================================================================
# PANEL MAKERS
# =================================================================

# 2x2 bar grid: BP-UP | BP-DOWN / KEGG-UP | KEGG-DOWN for a single contrast
make_quad_bars <- function(enr, title_str, out_stem, n_up, n_dn) {
  make_panel <- function(d, title_p, fill_col) {
    if (is.null(d) || !nrow(d)) {
      p <- ggplot() + annotate("text", x = 1, y = 1,
                               label = "No significant terms",
                               family = "Arial", size = 3,
                               color = "grey35") +
        theme_void(base_family = "Arial") +
        labs(title = title_p) +
        theme(plot.title = element_text(family = "Arial", face = "bold",
                                        size = 10))
      return(p)
    }
    d2 <- d %>% filter(!is.na(Description), Description != "", Description != "NA") %>%
      distinct(Description, .keep_all = TRUE) %>%
      head(N_TOP) %>%
      mutate(Description = substr(Description, 1, 45),
             Description = factor(Description, levels = rev(Description)),
             neglog = -log10(pmax(p.adjust, 1e-30)))
    ggplot(d2, aes(x = neglog, y = Description, fill = neglog)) +
      geom_col(alpha = 0.9, width = 0.72, color = "white") +
      geom_text(aes(label = sprintf("n=%d", Count)),
                hjust = -0.12, family = "Arial", size = 2.2,
                color = "grey15") +
      scale_fill_gradient(low = "#fef9e7", high = fill_col, guide = "none") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
      labs(x = expression(-log[10]*"(padj)"), y = NULL, title = title_p) +
      theme_grbv() +
      theme(axis.text.y = element_text(size = 6.3),
            plot.title  = element_text(size = 10, face = "bold"))
  }

  get_panel <- function(dir, onto) {
    if (is.null(enr)) return(NULL)
    enr %>% filter(direction == dir, ontology == onto) %>%
      arrange(p.adjust)
  }

  bpU  <- make_panel(get_panel("UP",   "GO:BP"),
                     sprintf("A  GO:BP  \u2014 UP (n = %d)",   n_up),
                     unname(GRBV_COLORS["up"]))
  bpD  <- make_panel(get_panel("DOWN", "GO:BP"),
                     sprintf("B  GO:BP  \u2014 DOWN (n = %d)", n_dn),
                     unname(GRBV_COLORS["down"]))
  kgU  <- make_panel(get_panel("UP",   "KEGG"),
                     "C  KEGG  \u2014 UP",
                     unname(GRBV_COLORS["up"]))
  kgD  <- make_panel(get_panel("DOWN", "KEGG"),
                     "D  KEGG  \u2014 DOWN",
                     unname(GRBV_COLORS["down"]))

  layout <- (bpU | bpD) / (kgU | kgD) +
    plot_layout(widths = c(1, 1), heights = c(1.25, 1)) +
    plot_annotation(
      title = title_str,
      subtitle = sprintf("Hypergeometric ORA, BH-FDR < 0.05;  UP: %d DEGs,  DOWN: %d DEGs;  top %d terms per panel",
                         n_up, n_dn, N_TOP),
      theme = theme(plot.title = element_text(family = "Arial", face = "bold",
                                              size = 12, hjust = 0.5,
                                              margin = margin(b = 2)),
                    plot.subtitle = element_text(family = "Arial", size = 8.5,
                                                 hjust = 0.5, color = "grey25",
                                                 margin = margin(b = 6),
                                                 lineheight = 1.1))
    )

  save_grbv(layout, file.path(OUT_DIR, out_stem),
            width_mm = GRBV$W_DOUBLE_MM, height_mm = 200)
}

# Stratified heatmap: pathways x (stratum, direction) for Q1 / Q2
make_stratum_heatmap <- function(per_stratum_de, stratum_order,
                                 title_str, subtitle_str, out_stem,
                                 ontology = "GO:BP") {
  # per_stratum_de: named list of DE data frames
  all_rows <- list()
  counts   <- list()
  for (nm in names(per_stratum_de)) {
    d <- per_stratum_de[[nm]]
    sets <- split_dir(d)
    counts[[nm]] <- list(UP = length(sets$UP), DOWN = length(sets$DOWN))
    for (dir in c("UP","DOWN")) {
      t2g <- if (ontology == "GO:BP") term2gene_bp else term2gene_kg
      t2n <- if (ontology == "GO:BP") term2name_bp else term2name_kg
      r <- enrich_one(sets[[dir]], t2g, t2n)
      if (!is.null(r) && nrow(r))
        all_rows[[length(all_rows)+1]] <- r %>%
          mutate(stratum = nm, direction = dir, ontology = ontology)
    }
  }
  if (length(all_rows) == 0) {
    message("[", out_stem, "] no enrichment found; skipping")
    return(invisible(NULL))
  }
  enr <- bind_rows(all_rows) %>%
    filter(!is.na(Description), Description != "", Description != "NA") %>%
    distinct(stratum, direction, Description, .keep_all = TRUE)

  top_terms <- enr %>%
    group_by(Description) %>%
    summarise(n_cells  = n(),
              min_padj = min(p.adjust),
              .groups  = "drop") %>%
    arrange(desc(n_cells), min_padj) %>%
    head(N_MASTER) %>% pull(Description)

  # Truncate Description consistently in BOTH enr and top_terms so
  # factor() levels match after substr. Length-55 is enough for every
  # GO:BP term in our universe but we apply it once here for safety.
  enr <- enr %>%
    filter(!is.na(Description), nzchar(Description), Description != "NA") %>%
    mutate(Description = substr(Description, 1, 55))

  top_terms <- substr(top_terms, 1, 55)
  top_terms <- top_terms[!is.na(top_terms) & top_terms != "NA" & nzchar(top_terms)]

  hmap <- enr %>%
    filter(Description %in% top_terms) %>%
    mutate(stratum     = factor(stratum,   levels = stratum_order),
           direction   = factor(direction, levels = c("UP","DOWN")),
           Description = factor(Description,
                                levels = rev(intersect(top_terms,
                                                       unique(Description)))),
           neglog = -log10(pmax(p.adjust, 1e-30))) %>%
    droplevels()

  max_neg <- max(hmap$neglog, na.rm = TRUE)
  txt_col <- ifelse(hmap$neglog > max_neg * 0.6, "white", "grey15")

  p <- ggplot(hmap, aes(x = stratum, y = Description, fill = neglog)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = Count),
              family = "Arial", size = 2.0, color = txt_col) +
    facet_wrap(~ direction, nrow = 1, scales = "fixed") +
    scale_fill_gradient(low = "#fef9e7",
                        high = unname(GRBV_COLORS["RB_pos"]),
                        name = expression(-log[10]*"(padj)")) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    labs(x = NULL, y = NULL,
         title = title_str, subtitle = subtitle_str) +
    theme_grbv() +
    theme(plot.title    = element_text(size = 12, face = "bold",
                                       hjust = 0.5),
          plot.subtitle = element_text(size = 8.5, hjust = 0.5,
                                       color = "grey25",
                                       margin = margin(b = 6),
                                       lineheight = 1.1),
          axis.text.x   = element_text(angle = 30, hjust = 1, size = 7.5),
          axis.text.y   = element_text(size = 6.5),
          strip.text    = element_text(face = "bold", size = 10),
          legend.key.height = unit(18, "pt"))

  save_grbv(p, file.path(OUT_DIR, out_stem),
            width_mm = GRBV$W_DOUBLE_MM, height_mm = 200)
  invisible(list(enr = enr, counts = counts))
}

# =================================================================
# S8  Q1 — age trajectory over 5 infection-onset cohorts (Q1b)
# =================================================================
q1_cohorts <- c("Q1b_cohort2016_vs_neg",
                "Q1b_cohort2019_vs_neg",
                "Q1b_cohort2020_vs_neg",
                "Q1b_cohort2021_vs_neg",
                "Q1b_cohort2022_vs_neg")
q1_labels  <- c("2016","2019","2020","2021","2022")
q1_de_list <- setNames(lapply(q1_cohorts, de_of), q1_labels)

s8_enr <- make_stratum_heatmap(
  q1_de_list, stratum_order = q1_labels,
  title_str    = "S8. Q1 pathway trajectory across infection-onset cohorts",
  subtitle_str = "GO:BP terms, direction-resolved; cell = -log10(padj), label = gene count",
  out_stem     = "S08_Q1_pathway_enrichment",
  ontology     = "GO:BP"
)

# =================================================================
# S9  Q2 — stage-resolved enrichment over 4 maturity stages
# =================================================================
q2_stages <- c("Pre-véraison","Véraison","Post-véraison","Harvest")
q2_de_list <- list(
  `Pre-véraison`  = de_of("Q2_Pre_veraison"),
  `Véraison`      = de_of("Q2_Veraison"),
  `Post-véraison` = de_of("Q2_Post_veraison"),
  `Harvest`       = de_of("Q2_Harvest")
)

s9_enr <- make_stratum_heatmap(
  q2_de_list, stratum_order = q2_stages,
  title_str    = "S9. Q2 pathway trajectory across maturity stages",
  subtitle_str = "GO:BP terms, direction-resolved; cell = -log10(padj), label = gene count",
  out_stem     = "S09_Q2_pathway_enrichment",
  ontology     = "GO:BP"
)

# =================================================================
# S10 Q3 — 2022 vs 2021 in RB+  (single contrast, UP/DOWN split)
# =================================================================
q3_de <- de_of("Q3_seasonRBpos")
q3_sets <- split_dir(q3_de)
make_quad_bars(
  enrich_by_direction(q3_de),
  "S10. Q3 pathway enrichment  \u2014  2022 vs 2021 in RB+",
  "S10_Q3_pathway_enrichment",
  length(q3_sets$UP), length(q3_sets$DOWN)
)

# =================================================================
# S11 Q4 — Season x Infection interaction (direction-preserving Q4d)
# =================================================================
q4_de <- de_of("Q4d_seasonXinfection_dir")
q4_sets <- split_dir(q4_de)
make_quad_bars(
  enrich_by_direction(q4_de),
  "S11. Q4 pathway enrichment  \u2014  Season \u00d7 Infection interaction",
  "S11_Q4_pathway_enrichment",
  length(q4_sets$UP), length(q4_sets$DOWN)
)

# =================================================================
# S12 Q5 — 2022 vs 2021 in RB-  (single contrast, UP/DOWN split)
# =================================================================
q5_de <- de_of("Q5_seasonRBneg")
q5_sets <- split_dir(q5_de)
make_quad_bars(
  enrich_by_direction(q5_de),
  "S12. Q5 pathway enrichment  \u2014  2022 vs 2021 in RB\u2212",
  "S12_Q5_pathway_enrichment",
  length(q5_sets$UP), length(q5_sets$DOWN)
)

# =================================================================
# S13 Q6 — pooled RB+ vs RB-  (single contrast, UP/DOWN split)
# =================================================================
q6_de <- de_of("Q6_pooled")
q6_sets <- split_dir(q6_de)
make_quad_bars(
  enrich_by_direction(q6_de),
  "S13. Q6 pathway enrichment  \u2014  pooled RB+ vs RB\u2212",
  "S13_Q6_pathway_enrichment",
  length(q6_sets$UP), length(q6_sets$DOWN)
)

# =================================================================
# S14 master landscape  Q1-Q6 x UP/DOWN
# =================================================================
rollup_de <- list(
  Q1 = do.call(rbind, q1_de_list) %>%
         group_by(gene) %>%
         slice_min(padj, n = 1, with_ties = FALSE) %>% ungroup() %>% as.data.frame(),
  Q2 = do.call(rbind, q2_de_list) %>%
         group_by(gene) %>%
         slice_min(padj, n = 1, with_ties = FALSE) %>% ungroup() %>% as.data.frame(),
  Q3 = q3_de, Q4 = q4_de, Q5 = q5_de, Q6 = q6_de
)

master_rows <- list()
for (qn in names(rollup_de)) {
  e <- enrich_by_direction(rollup_de[[qn]])
  if (!is.null(e) && nrow(e))
    master_rows[[length(master_rows)+1]] <- e %>%
      filter(ontology == "GO:BP") %>%
      mutate(question = qn)
}

if (length(master_rows)) {
  master_enr <- bind_rows(master_rows) %>%
    filter(!is.na(Description), Description != "", Description != "NA") %>%
    distinct(question, direction, Description, .keep_all = TRUE)

  master_top <- master_enr %>%
    group_by(Description) %>%
    summarise(n_cells  = n_distinct(paste(question, direction)),
              min_padj = min(p.adjust),
              .groups  = "drop") %>%
    arrange(desc(n_cells), min_padj) %>%
    head(30)

  master_enr <- master_enr %>%
    filter(!is.na(Description), nzchar(Description), Description != "NA") %>%
    mutate(Description = substr(Description, 1, 55))
  master_top_desc <- substr(master_top$Description, 1, 55)
  master_top_desc <- master_top_desc[!is.na(master_top_desc) &
                                     master_top_desc != "NA" &
                                     nzchar(master_top_desc)]

  master <- master_enr %>%
    filter(Description %in% master_top_desc) %>%
    mutate(question    = factor(question,  levels = names(rollup_de)),
           direction   = factor(direction, levels = c("UP","DOWN")),
           Description = factor(Description,
                                levels = rev(intersect(master_top_desc,
                                                       unique(Description)))),
           neglog = -log10(pmax(p.adjust, 1e-30))) %>%
    droplevels()

  p14 <- ggplot(master, aes(x = question, y = Description, fill = neglog)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = Count),
              family = "Arial", size = 2.1,
              color = ifelse(master$neglog > max(master$neglog) * 0.6,
                             "white","grey15")) +
    facet_wrap(~ direction, nrow = 1) +
    scale_fill_gradient(low = "#fef9e7",
                        high = unname(GRBV_COLORS["RB_pos"]),
                        name = expression(-log[10]*"(padj)")) +
    labs(x = "Research question", y = NULL,
         title = "S14. Master pathway landscape (top GO:BP terms)",
         subtitle = "30 most-recurrent terms across Q1-Q6, direction-resolved;  cell = -log10(padj), label = gene count") +
    theme_grbv() +
    theme(plot.title    = element_text(size = 12, face = "bold",
                                       hjust = 0.5),
          plot.subtitle = element_text(size = 8.5, hjust = 0.5,
                                       color = "grey25",
                                       margin = margin(b = 6)),
          axis.text.y   = element_text(size = 6.5),
          strip.text    = element_text(face = "bold", size = 10),
          legend.key.height = unit(20, "pt"))

  save_grbv(p14, file.path(OUT_DIR, "S14_master_pathway_landscape"),
            width_mm = GRBV$W_DOUBLE_MM, height_mm = 230)
} else {
  message("[supp_pathway] no enrichment found, skipping S14 master heatmap")
}

cat("[supp S8-S14] stratified + direction-resolved build complete.\n")
