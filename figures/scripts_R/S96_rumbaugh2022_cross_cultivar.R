# S96 — Cross-cultivar replication against Rumbaugh et al. 2022 Cabernet Sauvignon.
#
# Three panels:
#
#   A. DCL2 trajectory.
#      The endoribonuclease Dicer-like 2 (LOC100258190 / VIT_04s0023g00920)
#      log2-fold-change across the four ripening stages, separately for
#      Rumbaugh 2016 (cool), Rumbaugh 2017 (hot), our 2021 (cool), our 2022 (hot),
#      faceted by Rumbaugh's two rootstocks. Both papers show a veraison peak;
#      our data confirms it at padj 5.5e-11 (2021) / 2.3e-5 (2022).
#
#   B. Per-network concordance.
#      Spearman rank correlation between Rumbaugh's per-cell log2FC and our
#      Q2x per-stage-per-season log2FC, computed separately for each of the
#      three Rumbaugh co-expression networks and each of the two rootstocks.
#      Bar height = Spearman rho. Annotation = sign-agreement at our padj < 0.05.
#
#   C. Q1b cohort gradient.
#      Spearman rank correlation between Rumbaugh's harvest-mean per gene
#      (one value per rootstock) and each of our five infection-onset cohort
#      contrasts (cohort2016 = 5–6 years post-infection at our 2021 sampling,
#      down to cohort2022 = year 1). Peak agreement at cohort2020 places
#      Rumbaugh's vineyard at 1–2 years post-infection.
#
# Source CSVs (produced by farm_pipeline/scripts/rumbaugh2022_concordance.py):
#   rumbaugh2022_dcl2_trajectory.csv
#   rumbaugh2022_concordance_summary.csv
#   rumbaugh2022_q1b_concordance.csv

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(stringr)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ROOT    <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR <- file.path(ROOT, "figures/supp")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S96_rumbaugh2022_cross_cultivar")

dcl2 <- read_csv(file.path(ROOT, "rumbaugh2022_dcl2_trajectory.csv"),
                 show_col_types = FALSE)
summ <- read_csv(file.path(ROOT, "rumbaugh2022_concordance_summary.csv"),
                 show_col_types = FALSE)
q1b  <- read_csv(file.path(ROOT, "rumbaugh2022_q1b_concordance.csv"),
                 show_col_types = FALSE)

# Stage ordering for plotting
STAGE_ORDER <- c("Pre_veraison", "Veraison", "Post_veraison", "Harvest")
STAGE_LABS  <- c("Pre-V", "Veraison", "Post-V", "Harvest")
names(STAGE_LABS) <- STAGE_ORDER

# Network friendly labels for legend/axis
NETWORK_LABS <- c(
  S1_CaBP_hub  = "S1: Ca-binding hub",
  S2_DCL2_gain = "S2: DCL2-gain",
  S3_DCL2_loss = "S3: DCL2-loss"
)

# ---------------------------------------------------------------------------
# Panel A — DCL2 trajectory
# ---------------------------------------------------------------------------
# Keep one network (S2; DCL2 is identical between S2 and S3, but S2 lists it
# as the centralized DCL2-gain hub). Reshape to a four-series line plot per
# rootstock facet.

dcl2A <- dcl2 %>%
  filter(network == "S2_DCL2_gain") %>%
  mutate(stage = factor(stage, levels = STAGE_ORDER))

# Rumbaugh series: one trajectory per (rootstock, rumbaugh_year)
rumb_long <- dcl2A %>%
  distinct(rootstock, rumbaugh_year, stage, log2fc_rumbaugh) %>%
  mutate(study = "Rumbaugh 2022 CS",
         year_label = ifelse(rumbaugh_year == "2016",
                             "Rumbaugh 2016 (cool)",
                             "Rumbaugh 2017 (hot)")) %>%
  rename(log2fc = log2fc_rumbaugh)

# Our series: same trajectory; one value per (our_year, stage), independent
# of Rumbaugh's rootstock, so we drop rootstock duplicates.
ours_long <- dcl2A %>%
  distinct(our_year, stage, log2fc_ours_Q2x, padj_ours_Q2x) %>%
  mutate(study = "Merlot (this study)",
         year_label = ifelse(our_year == "2021",
                             "Merlot 2021 (cool)",
                             "Merlot 2022 (hot)")) %>%
  rename(log2fc = log2fc_ours_Q2x)

# Replicate the Merlot series across both Rumbaugh facets for visual comparison.
ours_long_facet <- bind_rows(
  ours_long %>% mutate(rootstock = "110R"),
  ours_long %>% mutate(rootstock = "420A")
)

trajA <- bind_rows(
  rumb_long  %>% select(rootstock, stage, study, year_label, log2fc),
  ours_long_facet %>% select(rootstock, stage, study, year_label, log2fc)
)

# Filter NAs (Post-veraison 2017 missing in Rumbaugh)
trajA <- trajA %>% filter(!is.na(log2fc))

# Color and linetype: cool blue-green vs hot magenta; Rumbaugh dashed, ours solid.
year_pal <- c(
  "Rumbaugh 2016 (cool)" = "#0F9D87",
  "Rumbaugh 2017 (hot)"  = "#C45BAA",
  "Merlot 2021 (cool)"   = "#0F9D87",
  "Merlot 2022 (hot)"    = "#C45BAA"
)
year_lty <- c(
  "Rumbaugh 2016 (cool)" = "dashed",
  "Rumbaugh 2017 (hot)"  = "dashed",
  "Merlot 2021 (cool)"   = "solid",
  "Merlot 2022 (hot)"    = "solid"
)
year_shape <- c(
  "Rumbaugh 2016 (cool)" = 1,
  "Rumbaugh 2017 (hot)"  = 1,
  "Merlot 2021 (cool)"   = 16,
  "Merlot 2022 (hot)"    = 16
)

pA <- ggplot(trajA, aes(x = stage, y = log2fc,
                        group = interaction(year_label),
                        colour = year_label,
                        linetype = year_label,
                        shape = year_label)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.3, fill = "white") +
  scale_x_discrete(labels = STAGE_LABS) +
  scale_colour_manual(values = year_pal, name = NULL,
                      breaks = c("Rumbaugh 2016 (cool)", "Rumbaugh 2017 (hot)",
                                 "Merlot 2021 (cool)", "Merlot 2022 (hot)")) +
  scale_linetype_manual(values = year_lty, name = NULL) +
  scale_shape_manual(values = year_shape, name = NULL) +
  facet_wrap(~ rootstock, ncol = 2,
             labeller = labeller(rootstock = c("110R" = "Rootstock: 110R",
                                                "420A" = "Rootstock: 420A"))) +
  labs(title = "(A) DCL2 antiviral switch peaks at veraison in both cultivars",
       subtitle = "Endoribonuclease Dicer-like 2 (LOC100258190). Cool seasons in green, hot seasons in magenta. Rumbaugh dashed, our Merlot solid.",
       x = "Berry development stage", y = "log2 fold-change (RB+ vs RB-)") +
  theme_grbv(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right",
        strip.text = element_text(face = "bold", size = 9))

# ---------------------------------------------------------------------------
# Panel B — per-network concordance
# ---------------------------------------------------------------------------
panelB_df <- summ %>%
  mutate(network = NETWORK_LABS[network],
         network = factor(network,
                          levels = unname(NETWORK_LABS)),
         # Signficant-gene sign-agreement annotation, as percent of agreeing
         # genes among genes with padj < 0.05 in our DESeq2 output
         annot = sprintf("rho = %.2f\n%.0f%% sign agree\n(n_sig = %d)",
                         spearman_rho,
                         100 * sign_agreement_sig,
                         n_sig_padj))

# Bar palette: distinct hues per network
network_pal <- c(
  "S1: Ca-binding hub" = "#1B7A9D",
  "S2: DCL2-gain"      = "#E1812B",
  "S3: DCL2-loss"      = "#8D6E2E"
)

pB <- ggplot(panelB_df, aes(x = network, y = spearman_rho,
                            fill = network)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.68, colour = "white", linewidth = 0.4,
           aes(group = rootstock, alpha = rootstock)) +
  geom_text(aes(label = annot, group = rootstock),
            position = position_dodge(width = 0.8),
            vjust = -0.18, size = 2.7, lineheight = 0.85,
            colour = "grey15") +
  scale_fill_manual(values = network_pal, guide = "none") +
  scale_alpha_manual(values = c("110R" = 1.0, "420A" = 0.55),
                     name = "Rumbaugh\nrootstock") +
  scale_y_continuous(limits = c(0, 0.55), breaks = seq(0, 0.5, 0.1),
                     expand = expansion(mult = c(0, 0.10))) +
  labs(title = "(B) All three networks replicate positively across both rootstocks",
       subtitle = "Spearman rank correlation (theirs vs ours) on per-cell log2FC; annotation = sign-agreement among our padj < 0.05 genes",
       x = "Rumbaugh 2022 differential co-expression network",
       y = "Spearman rho (per-cell log2FC)") +
  theme_grbv(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")

# ---------------------------------------------------------------------------
# Panel C — Q1b cohort gradient
# ---------------------------------------------------------------------------
panelC_df <- q1b %>%
  mutate(network = NETWORK_LABS[network],
         network = factor(network, levels = unname(NETWORK_LABS)),
         cohort_year = as.integer(str_extract(cohort, "\\d{4}")),
         years_since_2021 = 2021 - cohort_year,
         years_label = case_when(
           cohort_year == 2016 ~ "2016\n(5-6 y)",
           cohort_year == 2019 ~ "2019\n(2-3 y)",
           cohort_year == 2020 ~ "2020\n(1-2 y)",
           cohort_year == 2021 ~ "2021\n(0-1 y)",
           cohort_year == 2022 ~ "2022\n(year 1)"
         ),
         years_label = factor(years_label,
                              levels = c("2016\n(5-6 y)",
                                         "2019\n(2-3 y)",
                                         "2020\n(1-2 y)",
                                         "2021\n(0-1 y)",
                                         "2022\n(year 1)")))

pC <- ggplot(panelC_df,
             aes(x = years_label, y = rho_harvest,
                 colour = network, linetype = rootstock,
                 group = interaction(network, rootstock))) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = network_pal, name = NULL) +
  scale_linetype_manual(values = c("110R" = "solid", "420A" = "dashed"),
                        name = "Rumbaugh\nrootstock") +
  scale_y_continuous(limits = c(0, 0.75), breaks = seq(0, 0.7, 0.1)) +
  labs(title = "(C) Cross-paper agreement peaks at our 1-2 year cohort",
       subtitle = "Spearman rho on Rumbaugh-harvest-mean vs each of our infection-onset cohort contrasts. cohort2020 = 1-2 y post-infection at 2021 sampling.",
       x = "Our infection-onset cohort (year onset detected / years post-infection at sampling)",
       y = "Spearman rho (Rumbaugh harvest-mean vs Q1b)") +
  theme_grbv(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")

# ---------------------------------------------------------------------------
# Assemble
# ---------------------------------------------------------------------------
fig <- (pA / pB / pC) +
  plot_layout(heights = c(1.0, 0.95, 1.0)) +
  plot_annotation(
    title    = "S96 — Cross-cultivar replication against Rumbaugh 2022 Cabernet Sauvignon",
    subtitle = paste0(
      "Three layers of evidence. Panel A: the DCL2 veraison peak (their headline mechanism) replicates in our Merlot ",
      "(padj 5.5e-11 in 2021, 2.3e-5 in 2022). Panel B: all three of their co-expression networks carry positive ",
      "directional concordance with our data across both rootstocks, with 69-81% sign agreement among genes our ",
      "analysis calls significant. Panel C: agreement peaks at cohort2020 (1-2 y post-infection), placing Rumbaugh's ",
      "snapshot at the young-to-mid stage of the disease, not the chronic phase our oldest cohorts capture."),
    theme = theme(plot.title    = element_text(size = 12, face = "bold"),
                  plot.subtitle = element_text(size = 9, colour = "grey25",
                                               margin = margin(b = 6)))
  )

ggsave(paste0(OUT_STEM, ".png"), fig,
       width = GRBV$W_DOUBLE_MM, height = 290, units = "mm", dpi = 600)
ggsave(paste0(OUT_STEM, ".svg"), fig,
       width = GRBV$W_DOUBLE_MM, height = 290, units = "mm")

cat(sprintf(
  "[S96] wrote PNG + SVG to %s\n  Panels: A=DCL2 trajectory, B=per-network concordance, C=Q1b cohort gradient\n",
  OUT_STEM))
