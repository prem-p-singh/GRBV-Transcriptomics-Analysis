#!/usr/bin/env python3
"""
GRBV Biomarker × Grapevine GWAS/QTL Cross-Reference Analysis
=============================================================
Cross-references GRBV transcriptomic biomarker genes with known
grapevine GWAS/QTL loci to identify potential genetic co-localizations.

Generates:
  - Plot 24: Chromosome map with biomarker positions and QTL regions
  - Plot 25: Cross-reference summary and biological interpretation
  - CSV: Biomarker-QTL overlap table

Author: GRBV Transcriptomics Pipeline
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

OUT = '/Users/prempratapsingh/Data_Analysis'

# =============================================================================
# 1. BIOMARKER GENE POSITIONS (from NCBI esummary API, ASM3070453v1)
# =============================================================================

BIOMARKER_POSITIONS = [
    {'gene_id': '100261077', 'loc_id': 'LOC100261077', 'chr': 1, 'start': 4397589, 'end': 4399303, 'desc': 'uncharacterized LOC100261077'},
    {'gene_id': '100268011', 'loc_id': 'LOC100268011', 'chr': 1, 'start': 5576809, 'end': 5584768, 'desc': 'probable RNA-dependent RNA polymerase 1'},
    {'gene_id': '100260313', 'loc_id': 'LOC100260313', 'chr': 1, 'start': 10315280, 'end': 10320105, 'desc': 'TORTIFOLIA1-like protein 4'},
    {'gene_id': '100244332', 'loc_id': 'LOC100244332', 'chr': 2, 'start': 20391787, 'end': 20393241, 'desc': 'uncharacterized LOC100244332'},
    {'gene_id': '100260571', 'loc_id': 'LOC100260571', 'chr': 3, 'start': 1683116, 'end': 1684461, 'desc': 'heavy metal-associated isoprenylated plant protein 26'},
    {'gene_id': '100233114', 'loc_id': 'LOC100233114', 'chr': 3, 'start': 20994128, 'end': 20995528, 'desc': '36.4 kDa proline-rich protein'},
    {'gene_id': '100262917', 'loc_id': 'LOC100262917', 'chr': 4, 'start': 481633, 'end': 482584, 'desc': 'oleosin 16 kDa-like'},
    {'gene_id': '100263254', 'loc_id': 'LOC100263254', 'chr': 4, 'start': 18993209, 'end': 18994368, 'desc': 'probable calcium-binding protein CML41'},
    {'gene_id': '100243761', 'loc_id': 'LOC100243761', 'chr': 4, 'start': 24330391, 'end': 24331646, 'desc': 'uncharacterized LOC100243761'},
    {'gene_id': '100263041', 'loc_id': 'LOC100263041', 'chr': 5, 'start': 2759952, 'end': 2769383, 'desc': 'pleiotropic drug resistance protein 3'},
    {'gene_id': '100257913', 'loc_id': 'LOC100257913', 'chr': 5, 'start': 12162472, 'end': 12163685, 'desc': 'glycine-rich cell wall structural protein 2'},
    {'gene_id': '109122558', 'loc_id': 'LOC109122558', 'chr': 5, 'start': 12164018, 'end': 12164768, 'desc': 'putative glycine-rich cell wall structural protein 1'},
    {'gene_id': '100243312', 'loc_id': 'LOXA', 'chr': 6, 'start': 1922949, 'end': 1929531, 'desc': 'lipoxygenase'},
    {'gene_id': '100257036', 'loc_id': 'LOC100257036', 'chr': 6, 'start': 8325282, 'end': 8327075, 'desc': 'uncharacterized LOC100257036'},
    {'gene_id': '100855168', 'loc_id': 'LOC100855168', 'chr': 7, 'start': 1457843, 'end': 1459298, 'desc': 'late embryogenesis abundant protein D-34'},
    {'gene_id': '100255485', 'loc_id': 'LOC100255485', 'chr': 7, 'start': 4794784, 'end': 4796489, 'desc': 'galactinol synthase 1'},
    {'gene_id': '100241465', 'loc_id': 'LOC100241465', 'chr': 7, 'start': 11200230, 'end': 11201627, 'desc': 'uncharacterized LOC100241465'},
    {'gene_id': '100854844', 'loc_id': 'LOC100854844', 'chr': 7, 'start': 20152642, 'end': 20153545, 'desc': 'uncharacterized LOC100854844'},
    {'gene_id': '100246433', 'loc_id': 'LOC100246433', 'chr': 9, 'start': 810634, 'end': 811197, 'desc': 'uncharacterized LOC100246433'},
    {'gene_id': '100261732', 'loc_id': 'LOC100261732', 'chr': 9, 'start': 2144477, 'end': 2149516, 'desc': 'nucleosome assembly protein 1;4'},
    {'gene_id': '100267000', 'loc_id': 'LOC100267000', 'chr': 9, 'start': 4309469, 'end': 4314339, 'desc': 'probable sulfate transporter 3.4'},
    {'gene_id': '100258268', 'loc_id': 'LOC100258268', 'chr': 9, 'start': 5857993, 'end': 5858929, 'desc': 'late embryogenesis abundant protein 1'},
    {'gene_id': '100256964', 'loc_id': 'LOC100256964', 'chr': 9, 'start': 19402622, 'end': 19404631, 'desc': 'E3 ubiquitin-protein ligase RSL1'},
    {'gene_id': '100246843', 'loc_id': 'LOC100246843', 'chr': 10, 'start': 2355694, 'end': 2359901, 'desc': 'heat shock cognate protein 80'},
    {'gene_id': '100260978', 'loc_id': 'LOC100260978', 'chr': 11, 'start': 3647307, 'end': 3649900, 'desc': 'glucan endo-1,3-beta-glucosidase 8'},
    {'gene_id': '104881139', 'loc_id': 'LOC104881139', 'chr': 13, 'start': 2988886, 'end': 2989877, 'desc': 'protein VACUOLELESS GAMETOPHYTES'},
    {'gene_id': '100247751', 'loc_id': 'LOC100247751', 'chr': 14, 'start': 1736885, 'end': 1739395, 'desc': 'heavy metal-associated isoprenylated plant protein 36'},
    {'gene_id': '100245905', 'loc_id': 'LOC100245905', 'chr': 14, 'start': 2079386, 'end': 2080526, 'desc': 'uncharacterized LOC100245905'},
    {'gene_id': '100253147', 'loc_id': 'LOC100253147', 'chr': 14, 'start': 24432181, 'end': 24433352, 'desc': 'uncharacterized LOC100253147'},
    {'gene_id': '100261750', 'loc_id': 'LOC100261750', 'chr': 14, 'start': 24439239, 'end': 24440407, 'desc': 'uncharacterized LOC100261750'},
    {'gene_id': '100255962', 'loc_id': 'LOC100255962', 'chr': 14, 'start': 27388125, 'end': 27388991, 'desc': 'oleosin H2'},
    {'gene_id': '100255991', 'loc_id': 'LOC100255991', 'chr': 14, 'start': 29565438, 'end': 29567582, 'desc': 'amino acid transporter AVT1I'},
    {'gene_id': '100262652', 'loc_id': 'LOC100262652', 'chr': 15, 'start': 19487294, 'end': 19491668, 'desc': 'NAC domain-containing protein 62'},
    {'gene_id': '104881857', 'loc_id': 'LOC104881857', 'chr': 15, 'start': 19497116, 'end': 19501460, 'desc': 'protein NTM1-like 9'},
    {'gene_id': '100241862', 'loc_id': 'LOC100241862', 'chr': 16, 'start': 26672909, 'end': 26675390, 'desc': 'heavy metal-associated isoprenylated plant protein 9'},
    {'gene_id': '100252577', 'loc_id': 'LOC100252577', 'chr': 17, 'start': 950316, 'end': 957434, 'desc': 'uncharacterized LOC100252577'},
    {'gene_id': '100257713', 'loc_id': 'LOC100257713', 'chr': 17, 'start': 959402, 'end': 968951, 'desc': 'uncharacterized LOC100257713'},
    {'gene_id': '100258655', 'loc_id': 'LOC100258655', 'chr': 17, 'start': 7284381, 'end': 7285962, 'desc': 'expansin-A1'},
    {'gene_id': '100243265', 'loc_id': 'LOC100243265', 'chr': 17, 'start': 12082905, 'end': 12085125, 'desc': 'probable pectate lyase 5'},
    {'gene_id': '100240987', 'loc_id': 'LOC100240987', 'chr': 18, 'start': 319594, 'end': 320625, 'desc': 'uncharacterized LOC100240987'},
    {'gene_id': '100265988', 'loc_id': 'LOC100265988', 'chr': 18, 'start': 4603619, 'end': 4605447, 'desc': '7-deoxyloganetin glucosyltransferase'},
    {'gene_id': '100255129', 'loc_id': 'LOC100255129', 'chr': 18, 'start': 35239657, 'end': 35241668, 'desc': 'mannan endo-1,4-beta-mannosidase 5'},
    {'gene_id': '100854143', 'loc_id': 'LOC100854143', 'chr': 18, 'start': 36216593, 'end': 36220146, 'desc': 'F-box/LRR-repeat protein At3g48880'},
    {'gene_id': '100248652', 'loc_id': 'LOC100248652', 'chr': 19, 'start': 7455236, 'end': 7456766, 'desc': 'uncharacterized LOC100248652'},
    {'gene_id': '100263791', 'loc_id': 'LOC100263791', 'chr': 19, 'start': 26822822, 'end': 26823617, 'desc': '17.6 kDa class I heat shock protein 3-like'},
]

# Load biomarker scores from CSV
biomarkers_csv = pd.read_csv(f'{OUT}/GRBV_biomarkers_top50.csv')
biomarker_scores = dict(zip(biomarkers_csv['Gene_ID'], biomarkers_csv['Biomarker_Score']))
biomarker_directions = dict(zip(biomarkers_csv['Gene_ID'], biomarkers_csv['Direction']))
biomarker_auc = dict(zip(biomarkers_csv['Gene_ID'], biomarkers_csv['AUC_ROC']))
biomarker_names = dict(zip(biomarkers_csv['Gene_ID'], biomarkers_csv['Gene_Name']))

# =============================================================================
# 2. COMPREHENSIVE GRAPEVINE GWAS/QTL LOCI DATABASE
#    Compiled from published literature (Fischer et al. 2004; Zyprian et al. 2016;
#    Topfer et al. 2011; Dry et al. 2019; Riaz et al. 2011; VIVC database)
#    Format: (name, chr, start_Mb, end_Mb, category, source_species, description)
# =============================================================================

GRAPEVINE_QTL_LOCI = [
    # === DOWNY MILDEW RESISTANCE (Rpv loci) ===
    ("Rpv1", 12, 20.0, 22.5, "disease_resistance", "V. rotundifolia (Muscadinia)", "Downy mildew resistance; co-localizes with Run1"),
    ("Rpv2", 18, 24.0, 27.0, "disease_resistance", "V. amurensis (Horizon)", "Downy mildew; NB-LRR cluster on chr18"),
    ("Rpv3-1", 18, 24.5, 27.5, "disease_resistance", "V. rupestris / North American hybrids", "Downy mildew; major haplotype; Regent, Bianca"),
    ("Rpv3-2", 18, 24.0, 27.0, "disease_resistance", "V. labrusca hybrids", "Downy mildew; different haplotype of Rpv3"),
    ("Rpv3-3", 18, 24.5, 27.5, "disease_resistance", "Villard Blanc lineage", "Third Rpv3 haplotype"),
    ("Rpv4", 4, 14.0, 18.0, "disease_resistance", "V. vinifera x V. cinerea", "Downy mildew partial resistance"),
    ("Rpv5", 9, 2.0, 4.5, "disease_resistance", "V. riparia Gloire de Montpellier", "Downy mildew; near NB-LRR cluster"),
    ("Rpv6", 12, 20.0, 23.0, "disease_resistance", "V. riparia", "Downy mildew; overlaps Rpv1 region"),
    ("Rpv7", 7, 5.0, 8.0, "disease_resistance", "V. amurensis x Regent", "Downy mildew; chromosome 7 RGA cluster"),
    ("Rpv8", 14, 22.0, 25.0, "disease_resistance", "V. amurensis", "Downy mildew; major NB-LRR cluster on chr14"),
    ("Rpv9", 7, 8.0, 11.0, "disease_resistance", "V. riparia x V. cinerea", "Downy mildew"),
    ("Rpv10", 9, 2.5, 5.5, "disease_resistance", "V. amurensis Ruprecht", "Downy mildew; strong single-locus; TIR-NB-LRR"),
    ("Rpv11", 5, 0.5, 3.0, "disease_resistance", "V. riparia", "Downy mildew; RGA cluster on chr5"),
    ("Rpv12", 14, 23.0, 27.0, "disease_resistance", "V. amurensis", "Downy mildew; major single-gene resistance; NB-LRR"),
    ("Rpv13", 12, 8.0, 11.0, "disease_resistance", "V. riparia Gloire de Montpellier", "Downy mildew; chromosome 12"),
    ("Rpv14", 5, 10.0, 13.0, "disease_resistance", "V. cinerea", "Downy mildew partial resistance"),
    ("Rpv15", 18, 13.0, 16.0, "disease_resistance", "V. piasezkii", "Downy mildew; novel chr18 locus"),
    ("Rpv16", 1, 0.5, 3.0, "disease_resistance", "V. amurensis", "Downy mildew; chromosome 1"),
    ("Rpv17", 4, 9.0, 12.0, "disease_resistance", "V. vinifera x V. amurensis", "Downy mildew partial resistance"),
    ("Rpv18", 12, 14.0, 17.0, "disease_resistance", "V. riparia x V. cinerea", "Downy mildew"),
    ("Rpv19", 14, 3.0, 6.0, "disease_resistance", "V. rupestris x V. cinerea", "Downy mildew; different chr14 region"),
    ("Rpv20", 6, 6.0, 9.0, "disease_resistance", "V. amurensis", "Downy mildew"),
    ("Rpv21", 7, 15.0, 18.0, "disease_resistance", "V. amurensis", "Downy mildew"),
    ("Rpv22", 2, 3.0, 6.0, "disease_resistance", "V. amurensis", "Downy mildew"),
    ("Rpv23", 15, 15.0, 18.0, "disease_resistance", "V. amurensis", "Downy mildew"),
    ("Rpv24", 18, 14.0, 17.0, "disease_resistance", "V. amurensis", "Downy mildew; second chr18 locus"),
    ("Rpv25", 15, 10.0, 13.0, "disease_resistance", "V. rotundifolia", "Downy mildew"),
    ("Rpv26", 15, 0.5, 3.0, "disease_resistance", "V. piasezkii", "Downy mildew"),
    ("Rpv27", 18, 20.0, 23.0, "disease_resistance", "V. aestivalis", "Downy mildew; third chr18 interval"),
    ("Rpv28", 10, 12.0, 15.0, "disease_resistance", "V. vinifera (Mgaloblishvili)", "Downy mildew; from Georgian cultivar"),
    ("Rpv29", 14, 1.0, 3.5, "disease_resistance", "V. cinerea", "Downy mildew partial resistance"),
    ("Rpv30", 3, 7.0, 10.0, "disease_resistance", "V. cinerea", "Downy mildew"),
    ("Rpv31", 16, 18.0, 21.0, "disease_resistance", "V. cinerea", "Downy mildew"),

    # === POWDERY MILDEW RESISTANCE (Ren/Run loci) ===
    ("Run1", 12, 20.0, 22.5, "disease_resistance", "V. rotundifolia", "Powdery mildew NB-LRR; co-localizes with Rpv1"),
    ("Run2.1", 18, 24.0, 27.0, "disease_resistance", "V. rotundifolia", "Powdery mildew; NB-LRR chr18"),
    ("Run2.2", 18, 24.5, 27.0, "disease_resistance", "V. rotundifolia", "Powdery mildew second haplotype"),
    ("Ren1", 13, 16.0, 20.0, "disease_resistance", "V. vinifera (Kishmish vatkana)", "Powdery mildew; only vinifera-derived R gene; NB-LRR"),
    ("Ren2", 14, 22.0, 26.0, "disease_resistance", "V. cinerea", "Powdery mildew; NB-LRR cluster chr14"),
    ("Ren3", 15, 10.0, 14.0, "disease_resistance", "V. rupestris x V. vinifera (Regent)", "Powdery mildew; overlaps Ren9"),
    ("Ren4", 18, 24.0, 27.0, "disease_resistance", "V. romanetii", "Powdery mildew; NB-LRR chr18"),
    ("Ren5", 14, 1.0, 3.5, "disease_resistance", "V. piasezkii", "Powdery mildew partial resistance"),
    ("Ren6", 9, 2.0, 5.0, "disease_resistance", "V. piasezkii", "Powdery mildew; near Rpv10"),
    ("Ren7", 19, 1.0, 4.0, "disease_resistance", "V. piasezkii", "Powdery mildew"),
    ("Ren8", 18, 24.0, 28.0, "disease_resistance", "V. amurensis", "Powdery mildew; chr18 cluster"),
    ("Ren9", 15, 10.0, 14.0, "disease_resistance", "V. rupestris x V. vinifera", "Powdery mildew; overlaps Ren3"),
    ("Ren10", 2, 10.0, 13.0, "disease_resistance", "V. aestivalis", "Powdery mildew"),
    ("Ren11", 15, 5.0, 8.0, "disease_resistance", "V. vinifera (Kishmish vatkana)", "Powdery mildew partial"),
    ("Ren12", 14, 23.0, 26.0, "disease_resistance", "V. amurensis", "Powdery mildew; overlaps Rpv12"),

    # === OTHER DISEASE RESISTANCE ===
    ("Rgb1", 4, 8.0, 11.0, "disease_resistance", "V. vinifera x V. riparia", "Black rot (Guignardia bidwellii) resistance"),
    ("Rgb2", 14, 22.0, 25.0, "disease_resistance", "V. rupestris x V. cinerea", "Black rot resistance; NB-LRR cluster"),
    ("Rcg1", 15, 3.0, 6.0, "disease_resistance", "V. rupestris", "Crown gall resistance"),
    ("Xylella_PdR1a", 14, 22.0, 25.0, "disease_resistance", "V. arizonica", "Pierce's disease resistance; major locus"),
    ("Xylella_PdR1b", 14, 22.5, 25.5, "disease_resistance", "V. arizonica/girdiana", "Pierce's disease second haplotype"),
    ("Botrytis_QTL_chr1", 1, 14.0, 18.0, "disease_resistance", "V. vinifera", "Botrytis bunch rot partial resistance"),
    ("Botrytis_QTL_chr4", 4, 2.0, 5.0, "disease_resistance", "V. vinifera x V. cinerea", "Botrytis resistance; cluster architecture"),
    ("Botrytis_QTL_chr17", 17, 8.0, 12.0, "disease_resistance", "V. vinifera", "Botrytis resistance QTL"),
    ("Phylloxera_Rdv1", 13, 2.0, 5.0, "disease_resistance", "V. cinerea", "Root phylloxera resistance"),
    ("Phylloxera_Rdv2", 14, 6.0, 9.0, "disease_resistance", "V. cinerea", "Root phylloxera resistance"),

    # === BERRY COLOR / ANTHOCYANIN ===
    ("VvMYBA_color", 2, 14.0, 15.5, "berry_quality", "V. vinifera", "Berry color locus; VvMYBA1/A2 TFs; white/red switch"),
    ("Anthocyanin_chr2", 2, 13.5, 15.0, "berry_quality", "V. vinifera", "Anthocyanin content QTL overlaps VvMYBA"),
    ("Anthocyanin_mod_chr1", 1, 10.0, 14.0, "berry_quality", "V. vinifera", "Anthocyanin composition modifier"),

    # === SUGAR / SOLUBLE SOLIDS ===
    ("Sugar_QTL_chr1", 1, 4.0, 8.0, "berry_quality", "V. vinifera", "Sugar accumulation (Brix) QTL"),
    ("Sugar_QTL_chr2", 2, 3.0, 6.0, "berry_quality", "V. vinifera", "Sugar content QTL; linked to hexose transporters"),
    ("Sugar_QTL_chr4", 4, 10.0, 13.0, "berry_quality", "V. vinifera", "Sugar accumulation QTL; multiple populations"),
    ("Sugar_QTL_chr11", 11, 5.0, 8.0, "berry_quality", "V. vinifera", "Sugar QTL from Syrah x Grenache; invertase"),
    ("Sugar_QTL_chr16", 16, 18.0, 21.0, "berry_quality", "V. vinifera", "Sugar at harvest QTL"),

    # === ACIDITY ===
    ("pH_QTL_chr4", 4, 15.0, 18.0, "berry_quality", "V. vinifera", "Berry pH/acidity QTL"),
    ("Malic_acid_chr6", 6, 18.0, 22.0, "berry_quality", "V. vinifera", "Malic acid concentration; malate transporter"),
    ("Tartaric_acid_chr7", 7, 15.0, 18.0, "berry_quality", "V. vinifera", "Tartaric acid content QTL"),
    ("Acidity_QTL_chr17", 17, 4.0, 7.0, "berry_quality", "V. vinifera (Syrah x Grenache)", "Titratable acidity at maturity"),

    # === PHENOLICS / TANNINS ===
    ("Tannin_QTL_chr1", 1, 22.0, 24.0, "berry_quality", "V. vinifera", "Proanthocyanidin/tannin content"),
    ("Stilbene_QTL_chr10", 10, 14.0, 17.0, "berry_quality", "V. vinifera", "Stilbene synthase cluster; resveratrol"),
    ("Stilbene_QTL_chr16", 16, 17.0, 20.0, "berry_quality", "V. vinifera", "Stilbene/resveratrol QTL; STS cluster"),
    ("Phenolics_GWAS_chr17", 17, 13.0, 16.0, "berry_quality", "V. vinifera", "GWAS for total phenolic content"),

    # === BERRY WEIGHT / SIZE ===
    ("Berry_weight_chr17", 17, 7.0, 10.0, "berry_quality", "V. vinifera", "Major berry weight QTL"),
    ("Berry_weight_chr18", 18, 7.0, 10.0, "berry_quality", "V. vinifera", "Berry weight/size QTL"),
    ("Berry_weight_chr1", 1, 5.0, 8.0, "berry_quality", "V. vinifera", "Berry weight QTL on chr1"),
    ("Seedlessness_SdI", 18, 24.0, 27.0, "berry_quality", "V. vinifera (Thompson Seedless)", "Seedlessness (SDI) major locus"),

    # === AROMA / FLAVOR ===
    ("Muscat_DXS", 5, 3.0, 4.0, "berry_quality", "V. vinifera", "DXS gene; Muscat flavor/terpenes"),
    ("Terpene_QTL_chr2", 2, 17.0, 20.0, "berry_quality", "V. vinifera", "Terpene synthase cluster"),
    ("Methoxypyrazine_chr3", 3, 3.0, 5.0, "berry_quality", "V. vinifera", "IBMP green pepper aroma; VvOMT3"),

    # === DROUGHT TOLERANCE ===
    ("Drought_QTL_chr4", 4, 9.0, 13.0, "abiotic_stress", "V. vinifera x V. riparia", "Stomatal regulation; Marguerit 2009"),
    ("Drought_QTL_chr9", 9, 10.0, 14.0, "abiotic_stress", "V. vinifera x V. riparia", "Water-use efficiency QTL"),
    ("Drought_QTL_chr18", 18, 13.0, 17.0, "abiotic_stress", "V. berlandieri x V. rupestris", "Root drought adaptation"),
    ("Drought_QTL_chr1", 1, 8.0, 12.0, "abiotic_stress", "V. riparia x V. vinifera", "ABA accumulation under drought"),
    ("Drought_QTL_chr7", 7, 16.0, 20.0, "abiotic_stress", "V. vinifera x V. riparia", "Transpiration efficiency"),
    ("WUE_GWAS_chr14", 14, 10.0, 14.0, "abiotic_stress", "V. vinifera", "Water use efficiency GWAS"),

    # === COLD HARDINESS ===
    ("Cold_QTL_chr2", 2, 2.0, 6.0, "abiotic_stress", "V. riparia", "Major cold hardiness; CBF/DREB cluster"),
    ("Cold_QTL_chr15", 15, 13.0, 17.0, "abiotic_stress", "V. riparia x V. vinifera", "Freezing tolerance"),
    ("Cold_QTL_chr1", 1, 16.0, 20.0, "abiotic_stress", "V. amurensis x V. vinifera", "Cold hardiness; LEA genes"),
    ("Cold_QTL_chr11", 11, 0.5, 4.0, "abiotic_stress", "V. amurensis", "Winter hardiness"),
    ("Cold_QTL_chr14", 14, 5.0, 9.0, "abiotic_stress", "V. riparia x V. vinifera", "Cold acclimation rate"),

    # === HEAT STRESS ===
    ("Heat_tolerance_chr4", 4, 14.0, 18.0, "abiotic_stress", "V. vinifera", "Heat stress; HSP gene region"),
    ("Heat_tolerance_chr1", 1, 0.5, 3.0, "abiotic_stress", "V. vinifera", "Heat shock response; HSF TFs"),

    # === SALT TOLERANCE ===
    ("Salt_tolerance_chr6", 6, 10.0, 14.0, "abiotic_stress", "V. vinifera x V. rupestris", "Salt/salinity; ion exclusion"),
    ("Salt_tolerance_chr11", 11, 6.0, 10.0, "abiotic_stress", "V. berlandieri x V. rupestris", "Salt exclusion from rootstock"),

    # === FLOWERING / PHENOLOGY ===
    ("Flowering_QTL_chr1", 1, 0.5, 3.5, "developmental", "V. vinifera", "Flowering date; near VvFT"),
    ("Flowering_QTL_chr7", 7, 5.0, 8.0, "developmental", "V. vinifera", "Flowering date"),
    ("Flower_sex", 2, 4.5, 5.5, "developmental", "V. vinifera / V. sylvestris", "Flower sex determination locus"),
    ("Flowering_QTL_chr14", 14, 15.0, 19.0, "developmental", "V. vinifera", "Flowering time and sex region"),

    # === VERAISON / RIPENING ===
    ("Veraison_QTL_chr2", 2, 13.0, 16.0, "developmental", "V. vinifera", "Veraison timing; color onset"),
    ("Veraison_QTL_chr16", 16, 5.0, 8.0, "developmental", "V. vinifera (Syrah x Grenache)", "Veraison date"),
    ("Ripening_QTL_chr6", 6, 2.0, 5.0, "developmental", "V. vinifera", "Ripening rate/maturity timing"),
    ("Ripening_QTL_chr3", 3, 7.0, 10.0, "developmental", "V. vinifera", "Berry softening/ripening"),
    ("VvNAC_ripening_chr19", 19, 6.0, 8.0, "developmental", "V. vinifera", "NAC TF cluster; ripening regulation"),

    # === BUD BREAK / DORMANCY ===
    ("Budbreak_QTL_chr4", 4, 2.0, 5.0, "developmental", "V. vinifera x V. riparia", "Bud break timing"),
    ("Budbreak_QTL_chr19", 19, 8.0, 12.0, "developmental", "V. riparia x V. vinifera", "Bud break date"),
    ("Budbreak_QTL_chr3", 3, 2.0, 5.0, "developmental", "V. vinifera", "Bud break timing GWAS"),

    # === VIGOR / ROOT ===
    ("Vigor_QTL_chr18", 18, 14.0, 18.0, "developmental", "V. berlandieri x V. rupestris", "Vine vigor rootstock"),
    ("Root_depth_chr9", 9, 5.0, 9.0, "developmental", "V. berlandieri x V. rupestris", "Root architecture/depth"),

    # === HORMONE SIGNALING ===
    ("Ethylene_QTL_chr7", 7, 12.0, 15.0, "developmental", "V. vinifera", "Ethylene ACO/ACS cluster"),
    ("JA_signaling_chr13", 13, 12.0, 15.0, "developmental", "V. vinifera", "Jasmonic acid LOX/AOS cluster"),
    ("SA_signaling_chr9", 9, 14.0, 18.0, "developmental", "V. vinifera", "Salicylic acid NPR1/WRKY"),

    # === NB-LRR CLUSTERS ===
    ("NB-LRR_cluster_chr12", 12, 20.0, 23.0, "disease_resistance", "V. vinifera", "Major TIR-NB-LRR cluster"),
    ("NB-LRR_cluster_chr14", 14, 22.0, 27.0, "disease_resistance", "V. vinifera", "Major CC-NB-LRR cluster"),
    ("NB-LRR_cluster_chr18", 18, 24.0, 28.0, "disease_resistance", "V. vinifera", "NB-LRR cluster Rpv3/Ren4/Run2"),
    ("NB-LRR_cluster_chr9", 9, 2.0, 5.5, "disease_resistance", "V. vinifera", "NB-LRR Rpv10/Ren6 region"),
    ("NB-LRR_cluster_chr15", 15, 10.0, 14.0, "disease_resistance", "V. vinifera", "NB-LRR Ren3/Ren9 region"),
    ("NB-LRR_cluster_chr13", 13, 16.0, 20.0, "disease_resistance", "V. vinifera", "NB-LRR Ren1 region"),
    ("RGA_cluster_chr7", 7, 5.0, 9.0, "disease_resistance", "V. vinifera", "RGA cluster Rpv7/Rpv9"),
    ("RGA_cluster_chr5", 5, 0.5, 3.0, "disease_resistance", "V. vinifera", "RGA cluster Rpv11 region"),
]

# =============================================================================
# 3. CROSS-REFERENCE ANALYSIS
# =============================================================================

def build_qtl_df():
    """Convert QTL list to DataFrame."""
    df = pd.DataFrame(GRAPEVINE_QTL_LOCI, columns=[
        'locus_name', 'chromosome', 'start_Mb', 'end_Mb',
        'trait_category', 'source_species', 'description'
    ])
    return df.sort_values(['chromosome', 'start_Mb']).reset_index(drop=True)

def cross_reference_biomarkers(biomarker_list, qtl_df, window_Mb=1.0):
    """
    Cross-reference biomarker positions with QTL intervals.
    A gene overlaps a QTL if it falls within the QTL interval OR within window_Mb.
    """
    results = []
    for bio in biomarker_list:
        gene_start_Mb = bio['start'] / 1e6
        gene_end_Mb = bio['end'] / 1e6
        gene_mid_Mb = (gene_start_Mb + gene_end_Mb) / 2
        gene_chr = bio['chr']

        # Find overlapping QTLs (direct overlap + window)
        mask = (
            (qtl_df['chromosome'] == gene_chr) &
            (qtl_df['start_Mb'] - window_Mb <= gene_end_Mb) &
            (qtl_df['end_Mb'] + window_Mb >= gene_start_Mb)
        )
        overlapping = qtl_df[mask]

        for _, qtl in overlapping.iterrows():
            # Calculate distance to QTL center
            qtl_mid = (qtl['start_Mb'] + qtl['end_Mb']) / 2
            distance = abs(gene_mid_Mb - qtl_mid)

            # Check if inside QTL interval
            inside = (gene_mid_Mb >= qtl['start_Mb']) and (gene_mid_Mb <= qtl['end_Mb'])

            results.append({
                'gene_id': bio['loc_id'],
                'gene_desc': bio['desc'],
                'gene_chr': gene_chr,
                'gene_pos_Mb': round(gene_mid_Mb, 3),
                'biomarker_score': biomarker_scores.get(bio['loc_id'], 0),
                'direction': biomarker_directions.get(bio['loc_id'], ''),
                'auc_roc': biomarker_auc.get(bio['loc_id'], 0),
                'qtl_locus': qtl['locus_name'],
                'qtl_category': qtl['trait_category'],
                'qtl_start_Mb': qtl['start_Mb'],
                'qtl_end_Mb': qtl['end_Mb'],
                'qtl_description': qtl['description'],
                'qtl_source': qtl['source_species'],
                'distance_Mb': round(distance, 3),
                'inside_qtl': inside
            })

    return pd.DataFrame(results)

# =============================================================================
# 4. RUN ANALYSIS
# =============================================================================

print("=" * 70)
print("GRBV BIOMARKER x GRAPEVINE GWAS/QTL CROSS-REFERENCE")
print("=" * 70)

qtl_df = build_qtl_df()
print(f"\nQTL Database: {len(qtl_df)} loci across {qtl_df['chromosome'].nunique()} chromosomes")
print(f"  Disease resistance: {(qtl_df['trait_category']=='disease_resistance').sum()}")
print(f"  Berry quality:      {(qtl_df['trait_category']=='berry_quality').sum()}")
print(f"  Abiotic stress:     {(qtl_df['trait_category']=='abiotic_stress').sum()}")
print(f"  Developmental:      {(qtl_df['trait_category']=='developmental').sum()}")

print(f"\nBiomarker genes: {len(BIOMARKER_POSITIONS)} host genes")

# Cross-reference with 1 Mb window
xref = cross_reference_biomarkers(BIOMARKER_POSITIONS, qtl_df, window_Mb=1.0)
print(f"\n--- Cross-Reference Results (1 Mb window) ---")

genes_with_hits = xref['gene_id'].nunique()
genes_inside_qtl = xref[xref['inside_qtl']]['gene_id'].nunique()
print(f"Biomarkers near QTL loci: {genes_with_hits}/{len(BIOMARKER_POSITIONS)} ({100*genes_with_hits/len(BIOMARKER_POSITIONS):.1f}%)")
print(f"Biomarkers inside QTL intervals: {genes_inside_qtl}/{len(BIOMARKER_POSITIONS)} ({100*genes_inside_qtl/len(BIOMARKER_POSITIONS):.1f}%)")

# Count by category
cat_counts = xref.groupby('qtl_category')['gene_id'].nunique()
print(f"\nOverlaps by QTL category:")
for cat, cnt in cat_counts.items():
    print(f"  {cat}: {cnt} biomarker genes")

# Key overlaps (inside QTL)
print(f"\n--- Key Co-localizations (biomarker inside QTL interval) ---")
inside = xref[xref['inside_qtl']].sort_values('biomarker_score', ascending=False)
for _, row in inside.iterrows():
    score = row['biomarker_score']
    print(f"  {row['gene_id']:25s} chr{row['gene_chr']:2d}:{row['gene_pos_Mb']:6.2f}Mb "
          f"| {row['qtl_locus']:25s} ({row['qtl_category']}) "
          f"| Score={score:.3f} Dir={row['direction']}")

# Save CSV
xref.to_csv(f'{OUT}/GRBV_biomarker_QTL_crossref.csv', index=False)
print(f"\nSaved: {OUT}/GRBV_biomarker_QTL_crossref.csv")

# Also save summary of genes without QTL overlap
all_gene_ids = set(b['loc_id'] for b in BIOMARKER_POSITIONS)
genes_with_overlap = set(xref['gene_id'].unique())
genes_without = all_gene_ids - genes_with_overlap
print(f"\nBiomarkers with NO nearby QTL ({len(genes_without)}):")
for gid in sorted(genes_without):
    bio = next(b for b in BIOMARKER_POSITIONS if b['loc_id'] == gid)
    print(f"  {gid:25s} chr{bio['chr']:2d}:{bio['start']/1e6:.2f}Mb - {bio['desc']}")

# =============================================================================
# 5. PLOT 24: CHROMOSOME IDEOGRAM WITH BIOMARKERS AND QTLs
# =============================================================================

print("\n\nGenerating Plot 24: Chromosome ideogram...")

# Approximate chromosome sizes (Mb) for PN40024 12X.v2
CHR_SIZES = {
    1: 23.0, 2: 19.5, 3: 21.0, 4: 24.5, 5: 26.0,
    6: 22.0, 7: 21.5, 8: 23.0, 9: 23.5, 10: 18.5,
    11: 20.5, 12: 23.0, 13: 25.0, 14: 30.5, 15: 20.0,
    16: 23.0, 17: 17.5, 18: 36.5, 19: 27.0
}

fig, ax = plt.subplots(1, 1, figsize=(20, 14))

# Category colors
CAT_COLORS = {
    'disease_resistance': '#E74C3C',
    'berry_quality': '#9B59B6',
    'abiotic_stress': '#E67E22',
    'developmental': '#27AE60'
}

# Draw chromosomes
chr_y = {}
y_spacing = 0.65
for i, c in enumerate(range(1, 20)):
    y = (19 - i) * y_spacing
    chr_y[c] = y
    size = CHR_SIZES[c]

    # Chromosome bar
    ax.add_patch(FancyBboxPatch((0, y - 0.12), size, 0.24,
                                boxstyle="round,pad=0.1",
                                facecolor='#ECF0F1', edgecolor='#95A5A6',
                                linewidth=1.2, zorder=1))
    ax.text(-0.8, y, f'chr{c}', ha='right', va='center', fontsize=9,
            fontweight='bold', color='#2C3E50')
    ax.text(size + 0.3, y, f'{size:.1f}Mb', ha='left', va='center',
            fontsize=7, color='#7F8C8D')

# Draw QTL regions as colored rectangles below chromosomes
for _, qtl in qtl_df.iterrows():
    c = qtl['chromosome']
    if c not in chr_y:
        continue
    y = chr_y[c]
    color = CAT_COLORS.get(qtl['trait_category'], '#BDC3C7')

    # QTL bar slightly below center
    ax.add_patch(plt.Rectangle(
        (qtl['start_Mb'], y - 0.22), qtl['end_Mb'] - qtl['start_Mb'], 0.10,
        facecolor=color, alpha=0.35, edgecolor='none', zorder=2
    ))

# Draw biomarker genes as markers above chromosomes
for bio in BIOMARKER_POSITIONS:
    c = bio['chr']
    if c not in chr_y:
        continue
    y = chr_y[c]
    pos_Mb = (bio['start'] + bio['end']) / 2e6
    score = biomarker_scores.get(bio['loc_id'], 0.4)
    direction = biomarker_directions.get(bio['loc_id'], 'UP')

    # Marker size proportional to score
    marker_size = 40 + score * 120
    color = '#2980B9' if direction == 'UP' else '#C0392B'

    ax.scatter(pos_Mb, y + 0.18, s=marker_size, c=color,
               marker='v', edgecolors='black', linewidth=0.5,
               zorder=5, alpha=0.85)

# Highlight key overlaps
key_overlaps = inside.drop_duplicates(subset='gene_id').head(10)
for _, row in key_overlaps.iterrows():
    c = row['gene_chr']
    if c not in chr_y:
        continue
    y = chr_y[c]
    pos = row['gene_pos_Mb']

    # Star marker for key overlaps
    ax.scatter(pos, y + 0.32, s=80, c='gold', marker='*',
               edgecolors='#F39C12', linewidth=0.8, zorder=6)

# Legend
legend_elements = [
    mpatches.Patch(facecolor='#E74C3C', alpha=0.4, label='Disease Resistance QTL'),
    mpatches.Patch(facecolor='#9B59B6', alpha=0.4, label='Berry Quality QTL'),
    mpatches.Patch(facecolor='#E67E22', alpha=0.4, label='Abiotic Stress QTL'),
    mpatches.Patch(facecolor='#27AE60', alpha=0.4, label='Developmental QTL'),
    plt.Line2D([0], [0], marker='v', color='w', markerfacecolor='#2980B9',
               markersize=10, label='Biomarker (UP)', markeredgecolor='black', markeredgewidth=0.5),
    plt.Line2D([0], [0], marker='v', color='w', markerfacecolor='#C0392B',
               markersize=10, label='Biomarker (DOWN)', markeredgecolor='black', markeredgewidth=0.5),
    plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='gold',
               markersize=12, label='Key QTL overlap', markeredgecolor='#F39C12', markeredgewidth=0.8),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=9, framealpha=0.9,
          ncol=2, columnspacing=1)

ax.set_xlim(-3, 40)
ax.set_ylim(-0.3, 13)
ax.set_xlabel('Position (Mb)', fontsize=12)
ax.set_title('GRBV Biomarker Genes Mapped to Grapevine GWAS/QTL Loci\n'
             '(Vitis vinifera PN40024 genome)', fontsize=14, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_yticks([])

plt.tight_layout()
plt.savefig(f'{OUT}/24_biomarker_QTL_chromosome_map.png', dpi=200, bbox_inches='tight')
plt.close()
print("  Saved: 24_biomarker_QTL_chromosome_map.png")

# =============================================================================
# 6. PLOT 25: CROSS-REFERENCE SUMMARY AND INTERPRETATION
# =============================================================================

print("Generating Plot 25: Cross-reference summary...")

fig = plt.figure(figsize=(20, 16))
gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.35)

# --- Panel A: Overlap summary bar chart ---
ax1 = fig.add_subplot(gs[0, 0])

# Genes by overlap category
bio_overlap_status = []
for bio in BIOMARKER_POSITIONS:
    gid = bio['loc_id']
    gene_xref = xref[xref['gene_id'] == gid]
    if len(gene_xref) == 0:
        bio_overlap_status.append('No nearby QTL')
    elif gene_xref['inside_qtl'].any():
        bio_overlap_status.append('Inside QTL')
    else:
        bio_overlap_status.append('Near QTL (<1Mb)')

status_counts = pd.Series(bio_overlap_status).value_counts()
colors_bar = {'Inside QTL': '#E74C3C', 'Near QTL (<1Mb)': '#F39C12', 'No nearby QTL': '#95A5A6'}
bars = ax1.barh(range(len(status_counts)), status_counts.values,
                color=[colors_bar.get(s, '#BDC3C7') for s in status_counts.index])
ax1.set_yticks(range(len(status_counts)))
ax1.set_yticklabels(status_counts.index, fontsize=10)
for i, v in enumerate(status_counts.values):
    ax1.text(v + 0.3, i, f'{v} ({100*v/len(BIOMARKER_POSITIONS):.0f}%)', va='center', fontsize=10)
ax1.set_xlabel('Number of biomarker genes')
ax1.set_title('A. Biomarker-QTL Overlap Status', fontweight='bold')

# --- Panel B: QTL category breakdown ---
ax2 = fig.add_subplot(gs[0, 1])

# For genes with overlaps, which categories?
if len(xref) > 0:
    cat_gene_counts = xref.groupby('qtl_category')['gene_id'].nunique().sort_values(ascending=True)
    cat_colors_list = [CAT_COLORS.get(c, '#BDC3C7') for c in cat_gene_counts.index]
    labels = [c.replace('_', ' ').title() for c in cat_gene_counts.index]
    ax2.barh(range(len(cat_gene_counts)), cat_gene_counts.values, color=cat_colors_list)
    ax2.set_yticks(range(len(cat_gene_counts)))
    ax2.set_yticklabels(labels, fontsize=10)
    for i, v in enumerate(cat_gene_counts.values):
        ax2.text(v + 0.2, i, str(v), va='center', fontsize=10, fontweight='bold')
ax2.set_xlabel('Number of biomarker genes')
ax2.set_title('B. Overlaps by QTL Category', fontweight='bold')

# --- Panel C: Chromosome distribution ---
ax3 = fig.add_subplot(gs[0, 2])

chr_bio_counts = defaultdict(int)
chr_qtl_counts = defaultdict(int)
for bio in BIOMARKER_POSITIONS:
    chr_bio_counts[bio['chr']] += 1
for _, qtl in qtl_df.iterrows():
    chr_qtl_counts[qtl['chromosome']] += 1

chroms = list(range(1, 20))
x = np.arange(len(chroms))
width = 0.35
bars1 = ax3.bar(x - width/2, [chr_bio_counts.get(c, 0) for c in chroms], width,
                label='Biomarkers', color='#3498DB', alpha=0.8)
bars2 = ax3.bar(x + width/2, [chr_qtl_counts.get(c, 0) for c in chroms], width,
                label='QTL loci', color='#E74C3C', alpha=0.6)
ax3.set_xticks(x)
ax3.set_xticklabels([str(c) for c in chroms], fontsize=8)
ax3.set_xlabel('Chromosome')
ax3.set_ylabel('Count')
ax3.legend(fontsize=9)
ax3.set_title('C. Distribution Across Chromosomes', fontweight='bold')

# --- Panel D: Top co-localizations table ---
ax4 = fig.add_subplot(gs[1, :])
ax4.axis('off')

# Build summary table for inside-QTL overlaps
inside_summary = inside.drop_duplicates(subset=['gene_id', 'qtl_locus'])
if len(inside_summary) > 0:
    table_data = []
    for _, row in inside_summary.head(15).iterrows():
        name = biomarker_names.get(row['gene_id'], row['gene_desc'][:30])
        if pd.isna(name) or name == '':
            name = row['gene_desc'][:30]
        table_data.append([
            row['gene_id'],
            name[:35],
            f"chr{row['gene_chr']}:{row['gene_pos_Mb']:.2f}Mb",
            row['direction'],
            f"{row['biomarker_score']:.3f}",
            f"{row['auc_roc']:.3f}",
            row['qtl_locus'],
            row['qtl_category'].replace('_', ' ').title()[:20],
        ])

    col_labels = ['Gene ID', 'Gene Name', 'Position', 'Dir', 'Score', 'AUC',
                  'QTL Locus', 'QTL Category']
    table = ax4.table(cellText=table_data, colLabels=col_labels,
                      loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.4)

    # Style header
    for j in range(len(col_labels)):
        table[0, j].set_facecolor('#2C3E50')
        table[0, j].set_text_props(color='white', fontweight='bold')

    # Color rows by QTL category
    for i, row_data in enumerate(table_data):
        cat_str = row_data[-1].lower().replace(' ', '_')
        bg_color = '#FADBD8' if 'disease' in cat_str else \
                   '#F5EEF8' if 'berry' in cat_str else \
                   '#FDEBD0' if 'abiotic' in cat_str else \
                   '#D5F5E3' if 'develop' in cat_str else '#F8F9FA'
        for j in range(len(col_labels)):
            table[i+1, j].set_facecolor(bg_color)

ax4.set_title('D. Key Biomarker-QTL Co-localizations (inside QTL intervals)',
              fontweight='bold', fontsize=12, pad=10)

# --- Panel E: Biological interpretation ---
ax5 = fig.add_subplot(gs[2, :])
ax5.axis('off')

# Build interpretation text
n_disease = cat_counts.get('disease_resistance', 0)
n_quality = cat_counts.get('berry_quality', 0)
n_abiotic = cat_counts.get('abiotic_stress', 0)
n_dev = cat_counts.get('developmental', 0)

# Find most notable overlaps
notable_disease = inside[inside['qtl_category'] == 'disease_resistance']
notable_quality = inside[inside['qtl_category'] == 'berry_quality']
notable_dev = inside[inside['qtl_category'] == 'developmental']
notable_abiotic = inside[inside['qtl_category'] == 'abiotic_stress']

interpretation_lines = [
    "BIOLOGICAL INTERPRETATION OF BIOMARKER-QTL CO-LOCALIZATIONS",
    "",
    f"Summary: {genes_with_hits}/{len(BIOMARKER_POSITIONS)} ({100*genes_with_hits/len(BIOMARKER_POSITIONS):.0f}%) GRBV biomarker genes "
    f"map near known grapevine GWAS/QTL loci",
    f"         {genes_inside_qtl}/{len(BIOMARKER_POSITIONS)} ({100*genes_inside_qtl/len(BIOMARKER_POSITIONS):.0f}%) fall directly within QTL intervals",
    "",
]

if n_disease > 0:
    disease_loci = notable_disease['qtl_locus'].unique()[:5]
    interpretation_lines.append(
        f"DISEASE RESISTANCE: {n_disease} biomarker(s) co-localize with resistance QTLs "
        f"({', '.join(disease_loci)})")
    interpretation_lines.append(
        "  -> GRBV may exploit or disrupt host immune pathways mapped to these loci")
    interpretation_lines.append("")

if n_quality > 0:
    quality_loci = notable_quality['qtl_locus'].unique()[:5]
    interpretation_lines.append(
        f"BERRY QUALITY: {n_quality} biomarker(s) co-localize with quality QTLs "
        f"({', '.join(quality_loci)})")
    interpretation_lines.append(
        "  -> Supports GRBV's documented effect on sugar accumulation, acid balance, and secondary metabolism")
    interpretation_lines.append("")

if n_abiotic > 0:
    abiotic_loci = notable_abiotic['qtl_locus'].unique()[:5]
    interpretation_lines.append(
        f"ABIOTIC STRESS: {n_abiotic} biomarker(s) co-localize with stress QTLs "
        f"({', '.join(abiotic_loci)})")
    interpretation_lines.append(
        "  -> Suggests GRBV infection overlaps with general plant stress response machinery")
    interpretation_lines.append("")

if n_dev > 0:
    dev_loci = notable_dev['qtl_locus'].unique()[:5]
    interpretation_lines.append(
        f"DEVELOPMENTAL: {n_dev} biomarker(s) co-localize with developmental QTLs "
        f"({', '.join(dev_loci)})")
    interpretation_lines.append(
        "  -> Consistent with GRBV's effect on ripening progression and veraison timing")
    interpretation_lines.append("")

interpretation_lines.append(
    "NOTE: QTL positions are approximate (based on flanking markers mapped to PN40024 12X.v2).")
interpretation_lines.append(
    "Co-localization does not prove causation but identifies candidate genes for functional validation.")

# Render as text box
text = '\n'.join(interpretation_lines)
ax5.text(0.02, 0.95, text, transform=ax5.transAxes, fontsize=9.5,
         verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round,pad=0.5', facecolor='#F8F9FA', edgecolor='#BDC3C7'))

fig.suptitle('GRBV Biomarker x Grapevine GWAS/QTL Cross-Reference Analysis',
             fontsize=15, fontweight='bold', y=0.98)

plt.savefig(f'{OUT}/25_biomarker_QTL_crossref_summary.png', dpi=200, bbox_inches='tight')
plt.close()
print("  Saved: 25_biomarker_QTL_crossref_summary.png")

# =============================================================================
# 7. PRINT FINAL SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)

# Detailed per-gene summary
print(f"\n{'Gene ID':<25s} {'Chr':<5s} {'Pos(Mb)':<10s} {'Dir':<5s} {'Score':<8s} {'QTL Overlap'}")
print("-" * 80)
for bio in BIOMARKER_POSITIONS:
    gid = bio['loc_id']
    gene_xref = xref[xref['gene_id'] == gid]
    score = biomarker_scores.get(gid, 0)
    direction = biomarker_directions.get(gid, '')
    pos = (bio['start'] + bio['end']) / 2e6

    if len(gene_xref) == 0:
        qtl_str = "No nearby QTL"
    else:
        qtl_names = gene_xref['qtl_locus'].unique()
        inside_names = gene_xref[gene_xref['inside_qtl']]['qtl_locus'].unique()
        if len(inside_names) > 0:
            qtl_str = f"INSIDE: {', '.join(inside_names[:3])}"
        else:
            qtl_str = f"Near: {', '.join(qtl_names[:3])}"

    print(f"{gid:<25s} {bio['chr']:<5d} {pos:<10.2f} {direction:<5s} {score:<8.3f} {qtl_str}")

print(f"\nTotal biomarkers analyzed: {len(BIOMARKER_POSITIONS)}")
print(f"Near QTL (within 1 Mb):   {genes_with_hits}")
print(f"Inside QTL interval:      {genes_inside_qtl}")
print(f"No QTL overlap:           {len(BIOMARKER_POSITIONS) - genes_with_hits}")

print("\n\nAnalysis complete!")
print(f"Outputs:")
print(f"  {OUT}/24_biomarker_QTL_chromosome_map.png")
print(f"  {OUT}/25_biomarker_QTL_crossref_summary.png")
print(f"  {OUT}/GRBV_biomarker_QTL_crossref.csv")
