# Supplementary: mediation + chromosome + WGCNA detail — S49..S55.
#
#   S49  Mediation path diagrams for top 6 mediators (full detail)
#   S50  Direct vs indirect scatter + category bar (standalone)
#   S51  Chromosome DEG-density ideogram
#   S52  DEG enrichment Manhattan (full)
#   S53  WGCNA module-trait correlations (all 13 modules)
#   S54  Module eigengene profiles across stages
#   S55  Hub gene heatmaps for the top-4 infection-correlated modules

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table); library(ggrepel)
})
HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))
DATA_ROOT <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR   <- file.path(DATA_ROOT, "figures", "supp")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
stages <- GRBV_STAGES

# -----------------------------------------------------------------
# Mediation inputs
# -----------------------------------------------------------------
med <- fread(file.path(DATA_ROOT, "mediation_results.csv")) %>% as.data.frame()
ann <- fread(file.path(DATA_ROOT,
                       "farm_pipeline/results/fix11_annotations/gene_annotations.csv")) %>%
  as.data.frame()
name_map <- setNames(ann$product, ann$gene_id)

short <- c(
  LOC100249249 = "MBF1c",  LOC100853778 = "PSK6",
  LOC100261367 = "PP2-B15", LOC100253072 = "DREB1E",
  LOC100247528 = "GRX-C9",  LOC100232974 = "SPS1"
)

# -----------------------------------------------------------------
# S49 — top-6 mediator path diagrams
# -----------------------------------------------------------------
sig <- med %>% filter(sig_joint) %>% arrange(path_b_padj) %>% head(6)

mk_panel <- function(row) {
  gname <- name_map[row$gene]; if (is.na(gname)) gname <- ""
  gname <- substr(gname, 1, 38)
  g_short <- if (row$gene %in% names(short)) short[[row$gene]] else row$gene

  nodes <- data.frame(
    x = c(1.9, 5.0, 8.1), y = c(4.75, 8.15, 4.75),
    w = c(3.2, 2.8, 3.2), h = c(1.4, 1.4, 1.4),
    fill = c(unname(GRBV_COLORS["med_temp"]),
             unname(GRBV_COLORS["med_titer"]),
             unname(GRBV_COLORS["highlight"])),
    lab  = c("Temp", "Titer", g_short))
  arrs <- data.frame(
    x = c(3.4, 6.6, 3.5),   y = c(5.45, 7.75, 4.75),
    xend = c(4.2, 5.8, 6.5),yend = c(7.75, 5.45, 4.75),
    type = c("a","b","cprime"))

  ggplot() +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 10), expand = FALSE) +
    geom_rect(data = nodes, aes(xmin = x - w/2, xmax = x + w/2,
                                ymin = y - h/2, ymax = y + h/2),
              fill = nodes$fill, color = "black", linewidth = 0.4) +
    geom_text(data = nodes, aes(x = x, y = y, label = lab),
              family = "Arial", color = "white",
              fontface = "bold", size = 3.2) +
    geom_segment(data = arrs[arrs$type %in% c("a","b"), ],
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color = "black", linewidth = 0.5,
                 arrow = arrow(length = unit(6, "pt"), type = "closed")) +
    geom_segment(data = arrs[arrs$type == "cprime", ],
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color = unname(GRBV_COLORS["RB_pos"]), linewidth = 0.5,
                 linetype = "dashed",
                 arrow = arrow(length = unit(6, "pt"), type = "closed")) +
    annotate("text", x = 5, y = 3.3,
             label = sprintf("c = %+.3f   p = %s",
                             row$total_c, scales::scientific(row$total_p, 2)),
             family = "Arial", size = 2.7, hjust = 0.5) +
    annotate("text", x = 5, y = 2.5,
             label = sprintf("a x b = %+.3f   joint padj = %s",
                             row$indirect_ab,
                             scales::scientific(row$path_b_padj, 2)),
             family = "Arial", size = 2.7, hjust = 0.5) +
    annotate("text", x = 5, y = 1.7,
             label = sprintf("%% mediated = %+.1f%%",
                             row$prop_mediated * 100),
             family = "Arial", size = 3.0, hjust = 0.5,
             fontface = "bold", color = unname(GRBV_COLORS["highlight"])) +
    annotate("text", x = 5, y = 0.8, label = gname,
             family = "Arial", size = 2.4, hjust = 0.5,
             fontface = "italic", color = "grey25") +
    labs(title = row$gene) +
    theme_void(base_family = "Arial") +
    theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5,
                                    family = "Arial"))
}

panels <- lapply(seq_len(nrow(sig)), function(i) mk_panel(sig[i, ]))
S49 <- wrap_plots(panels, ncol = 3) +
  plot_annotation(
    title = "Supplementary S49. Mediation path diagrams for the top-6 joint-significant mediators",
    caption = "tv = 2,000 most variable genes; BH FDR<0.05 on Temp->Gene and Titer->Gene; |a x b| > 0.05",
    theme = theme(plot.title = element_text(family = "Arial", face = "bold",
                                            size = 12, hjust = 0.5),
                  plot.caption = element_text(family = "Arial", size = 8,
                                              color = "grey35", hjust = 0.5)))
save_grbv(S49, file.path(OUT_DIR, "S49_mediation_diagrams"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 140)

# -----------------------------------------------------------------
# S50 — direct vs indirect scatter + category bar
# -----------------------------------------------------------------
med50 <- med %>%
  mutate(category = case_when(
    sig_joint ~ "Titer-mediated",
    total_padj < 0.05 ~ "Temperature-only",
    TRUE ~ "Not significant"))
cat_color50 <- pal_mediation

lab50 <- med %>% filter(sig_joint) %>% arrange(path_b_padj) %>% head(8) %>%
  mutate(label = ifelse(gene %in% names(short), short[gene], gene))

pA50 <- ggplot(med50, aes(x = direct_cprime, y = indirect_ab, color = category)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_point(alpha = 0.55, size = 0.8) +
  geom_text_repel(data = lab50, aes(label = label),
                  color = "black", family = "Arial", fontface = "bold",
                  size = 2.5, box.padding = 0.4, segment.size = 0.2) +
  scale_color_manual(values = cat_color50, name = NULL) +
  labs(x = "Direct effect c'", y = "Indirect effect a x b",
       title = "A. Direct vs indirect temperature effects") +
  theme_grbv() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = c(0.02, 0.98), legend.justification = c(0, 1))

counts50 <- med50 %>% count(category) %>%
  mutate(category = factor(category,
                           levels = c("Titer-mediated","Temperature-only","Not significant")))
pB50 <- ggplot(counts50, aes(x = n, y = category, fill = category)) +
  geom_col(alpha = 0.92, width = 0.7, color = "white") +
  geom_text(aes(label = scales::comma(n)), hjust = -0.15,
            family = "Arial", size = 2.8) +
  scale_fill_manual(values = cat_color50, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20)),
                     labels = scales::label_comma()) +
  scale_y_discrete(limits = rev) +
  labs(x = "Number of genes", y = NULL,
       title = "B. Gene classification") +
  theme_grbv() +
  theme(plot.title = element_text(size = 10, face = "bold"))

S50 <- pA50 | pB50 +
  plot_annotation(
    title = "Supplementary S50. Direct vs titer-mediated temperature effects",
    theme = theme(plot.title = element_text(family = "Arial", face = "bold",
                                            size = 12, hjust = 0.5)))
save_grbv(S50, file.path(OUT_DIR, "S50_direct_vs_indirect_effects"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 110)

# -----------------------------------------------------------------
# S51, S52 — chromosome DEG ideogram + Manhattan
# -----------------------------------------------------------------
deg_dens <- fread(file.path(DATA_ROOT, "chromosome_deg_density.csv")) %>% as.data.frame()
deg_dens$chr <- factor(deg_dens$chr, levels = sort(unique(deg_dens$chr)))
qtl <- fread(file.path(DATA_ROOT, "GRBV_biomarker_QTL_crossref.csv")) %>% as.data.frame()

# S51 — ideogram (per-chromosome DEG density as horizontal bar, stacked by window)
dens_plot <- deg_dens %>%
  mutate(window_mid_mb = window_mid / 1e6)

S51 <- ggplot(dens_plot, aes(x = window_mid_mb, y = chr,
                             fill = deg_fraction)) +
  geom_tile(height = 0.8, width = 1.0, color = NA) +
  scale_fill_gradient(low = "#f5f5f5", high = unname(GRBV_COLORS["RB_pos"]),
                      name = "DEG fraction\nper 1 Mb window",
                      labels = scales::percent_format(accuracy = 1)) +
  scale_y_discrete(limits = rev) +
  labs(x = "Genomic position (Mb)", y = "Chromosome",
       title = "Supplementary S51. Chromosome DEG-density ideogram (1 Mb windows)") +
  theme_grbv() +
  theme(plot.title = element_text(size = 11, face = "bold"),
        legend.key.height = unit(18, "pt"))
save_grbv(S51, file.path(OUT_DIR, "S51_chromosome_deg_ideogram"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 130)

# S52 — Manhattan per-chromosome, stacked facet
# Cumulative position
chr_lens <- deg_dens %>% group_by(chr) %>%
  summarise(maxpos = max(window_end), .groups = "drop") %>%
  arrange(as.integer(as.character(chr))) %>%
  mutate(offset = lag(cumsum(maxpos), default = 0))
manhattan <- deg_dens %>% left_join(chr_lens %>% select(chr, offset), by = "chr") %>%
  mutate(cum_pos = window_mid + offset,
         neglog = -log10(pmax(enrichment_padj, 1e-30)))
chr_centers <- manhattan %>% group_by(chr) %>%
  summarise(center = (min(cum_pos) + max(cum_pos)) / 2, .groups = "drop")
# Alternating greys for chromosome strips (keeps red/blue reserved for DEG direction)
chr_cols <- rep(c("grey60", "grey85"),
                length.out = length(unique(manhattan$chr)))
names(chr_cols) <- levels(manhattan$chr)

S52 <- ggplot(manhattan, aes(x = cum_pos, y = neglog, color = chr)) +
  geom_point(alpha = 0.8, size = 0.8) +
  geom_hline(yintercept = -log10(0.05), color = "grey50",
             linetype = "dashed", linewidth = 0.3) +
  scale_color_manual(values = chr_cols, guide = "none") +
  scale_x_continuous(breaks = chr_centers$center,
                     labels = chr_centers$chr,
                     expand = c(0.01, 0.01)) +
  labs(x = "Chromosome",
       y = expression(-log[10]*"(window padj)"),
       title = "Supplementary S52. Per-chromosome DEG-enrichment Manhattan plot") +
  theme_grbv() +
  theme(plot.title = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 7),
        panel.grid.major.x = element_blank())
save_grbv(S52, file.path(OUT_DIR, "S52_deg_manhattan_plot"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 90)

# -----------------------------------------------------------------
# WGCNA
# -----------------------------------------------------------------
mt  <- fread(file.path(DATA_ROOT, "module_trait_correlations.csv")) %>% as.data.frame()
colnames(mt)[1] <- "module"
hub <- fread(file.path(DATA_ROOT, "hub_genes_per_module.csv")) %>% as.data.frame()

# S53 — Module-trait correlation heatmap (all 13 modules)
mt_long <- mt %>%
  pivot_longer(-module, names_to = "trait", values_to = "r") %>%
  mutate(trait = recode(trait, "InfectionYear" = "Years infected"))
ord <- mt %>% mutate(ai = abs(Infection)) %>% arrange(desc(ai)) %>% pull(module)
mt_long$module <- factor(mt_long$module, levels = rev(ord))
mt_long$trait  <- factor(mt_long$trait,
                         levels = c("Stage","Season","Infection","Years infected"))

S53 <- ggplot(mt_long, aes(x = trait, y = module, fill = r)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.2f", r)),
            family = "Arial", size = 2.5,
            color = ifelse(abs(mt_long$r) > 0.5, "white", "grey15")) +
  scale_fill_gradient2(low = unname(GRBV_COLORS["down"]),
                       mid = "white", high = unname(GRBV_COLORS["up"]),
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  labs(x = NULL, y = NULL,
       title = "Supplementary S53. WGCNA module-trait correlations (all 13 modules)") +
  theme_grbv() +
  theme(plot.title = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 8.5),
        axis.text.y = element_text(size = 8.5),
        legend.key.height = unit(20, "pt"))
save_grbv(S53, file.path(OUT_DIR, "S53_module_trait_correlations"),
          width_mm = GRBV$W_ONEHALF_MM, height_mm = 140)

# S54 — module eigengene profiles
MEs_f <- file.path(DATA_ROOT, "MEs.csv")
if (file.exists(MEs_f)) {
  MEs <- fread(MEs_f) %>% as.data.frame()
  # First col is sample; keep only ME* columns for the pivot
  names(MEs)[1] <- "Sample"
  me_cols <- grep("^ME", names(MEs), value = TRUE)
  MEs <- MEs[, c("Sample", me_cols)]
  meta <- read_tsv(file.path(DATA_ROOT, "sample_metadata.txt"),
                   show_col_types = FALSE) %>% as.data.frame()
  MEs_long <- MEs %>%
    pivot_longer(-Sample, names_to = "module", values_to = "eigen") %>%
    mutate(module = sub("^ME", "", module)) %>%
    inner_join(meta, by = "Sample") %>%
    mutate(Maturity_stage = factor(Maturity_stage, levels = stages),
           module = factor(module, levels = ord))
  summ <- MEs_long %>%
    group_by(module, Maturity_stage, Infection_Group) %>%
    summarise(mean = mean(eigen), se = sd(eigen)/sqrt(n()), .groups = "drop")

  S54 <- ggplot(summ, aes(x = Maturity_stage, y = mean,
                          color = Infection_Group,
                          group = Infection_Group)) +
    geom_line(linewidth = 0.6) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  width = 0.15, linewidth = 0.4) +
    geom_point(size = 1.4) +
    facet_wrap(~ module, ncol = 4) +
    scale_color_manual(values = GRBV_INFECTION_COLOR, name = NULL) +
    scale_x_discrete(labels = substr(GRBV_STAGE_LABEL, 1, 4)) +
    labs(x = NULL, y = "Module eigengene (mean)",
         title = "Supplementary S54. Module eigengene profiles across stages") +
    theme_grbv() +
    theme(plot.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5),
          strip.text  = element_text(size = 8, face = "bold"),
          legend.position = "top")
  save_grbv(S54, file.path(OUT_DIR, "S54_module_eigengene_profiles"),
            width_mm = GRBV$W_DOUBLE_MM, height_mm = 180)
}

# S55 — hub gene heatmaps for top-4 infection-correlated modules
top_mods <- ord[1:min(4, length(ord))]
hubs55 <- hub %>% filter(module %in% top_mods) %>%
  group_by(module) %>% arrange(desc(kME)) %>% slice_head(n = 5) %>%
  ungroup()

# Need lcpm for these genes
counts <- fread(file.path(DATA_ROOT, "counts_matrix.txt"),
                sep = "\t", header = TRUE, data.table = FALSE)
rownames(counts) <- counts[,1]; counts[,1] <- NULL
counts <- as.matrix(counts); storage.mode(counts) <- "integer"
keep <- rowSums(counts >= 10) >= 0.10 * ncol(counts); counts <- counts[keep, ]
lib <- colSums(counts); lcpm <- log2(sweep(counts, 2, lib, "/") * 1e6 + 1)
meta <- read_tsv(file.path(DATA_ROOT, "sample_metadata.txt"),
                 show_col_types = FALSE) %>% as.data.frame()
rownames(meta) <- meta$Sample; meta <- meta[colnames(lcpm), ]
meta$Maturity_stage <- factor(meta$Maturity_stage, levels = stages)

hub_genes <- intersect(unique(hubs55$gene), rownames(lcpm))
if (length(hub_genes) >= 4) {
  # Mean per stage x infection for heatmap
  Zmean <- lapply(hub_genes, function(g) {
    df <- data.frame(sample = colnames(lcpm), y = as.numeric(lcpm[g, ]))
    df <- merge(df, meta[, c("Sample","Maturity_stage","Infection_Group")],
                by.x = "sample", by.y = "Sample")
    df$group <- paste(df$Infection_Group, df$Maturity_stage, sep = " | ")
    aggregate(y ~ group, data = df, FUN = mean)
  })
  names(Zmean) <- hub_genes
  hm <- do.call(rbind, lapply(hub_genes, function(g) {
    d <- Zmean[[g]]; d$gene <- g; d
  }))
  # Z-score per gene across groups
  hm <- hm %>% group_by(gene) %>%
    mutate(z = (y - mean(y)) / sd(y)) %>% ungroup()

  # Attach module for faceting
  hm <- hm %>% left_join(hubs55 %>% select(gene, module), by = "gene") %>%
    mutate(module = factor(module, levels = top_mods),
           group  = factor(group, levels = c(paste("RB(-)", levels(meta$Maturity_stage)),
                                             paste("RB(+)", levels(meta$Maturity_stage)))))

  S55 <- ggplot(hm, aes(x = group, y = gene, fill = z)) +
    geom_tile(color = "white", linewidth = 0.3) +
    facet_wrap(~ module, ncol = 2, scales = "free_y") +
    scale_fill_gradient2(low = unname(GRBV_COLORS["down"]),
                         mid = "white", high = unname(GRBV_COLORS["up"]),
                         midpoint = 0, name = "z-score", limits = c(-2.5, 2.5),
                         oob = scales::squish) +
    labs(x = NULL, y = NULL,
         title = "Supplementary S55. Top-5 hub gene expression (z-score) for the top-4 infection-correlated modules") +
    theme_grbv() +
    theme(plot.title = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
          axis.text.y = element_text(size = 6.5),
          strip.text = element_text(size = 9, face = "bold"),
          legend.key.height = unit(18, "pt"))
  save_grbv(S55, file.path(OUT_DIR, "S55_hub_gene_heatmaps"),
            width_mm = GRBV$W_DOUBLE_MM, height_mm = 150)
}

cat("[supp_mediation_chr_wgcna] wrote S49..S55\n")
