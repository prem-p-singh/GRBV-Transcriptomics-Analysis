# `data/literature/` — published Vitis DEG tables for the comparative meta-analysis

Rules 11 and 12 of the FARM pipeline consume published DEG tables from this
directory. Each file is one comparison (one stress vs its own control) from
one paper. The pipeline maps heterogeneous identifiers to NCBI `LOC_` IDs,
computes overlap statistics against the GRBV DEG signature, and classifies
each GRBV DEG as GRBV-unique / viral-shared / pathogen-generic / stress-
generic.

---

## 1. Tier 1 study list (verified 2026-04-23)

Six papers ingested. Each citation verified on PubMed / DOI resolver; three
had to be corrected from earlier drafts after strict verification.

| Stress axis | Study | DOI / PMID | Platform, tissue | Notes |
|---|---|---|---|---|
| **GRBV** (same virus, different cultivar) | Blanco-Ulate B, Hopfer H, Figueroa-Balderas R, Ye Z, Rivero RM, Albacete A, Pérez-Alfocea F, Koyama R, Anderson MM, Smith RJ, Ebeler SE, Cantu D. 2017. Red blotch disease alters grape berry development and metabolism by interfering with the transcriptional and hormonal regulation of ripening. *J Exp Bot* 68(5):1225–1238. | 10.1093/jxb/erw506 / 28338755 | RNA-Seq; Zinfandel berries, two vineyards (Healdsburg, Oakville), four stages | 8 contrasts (2 vineyards × 4 stages). Full-panel table. |
| **GLRaV-3** (virus, different family) | Song YS, Hanner RH, Meng B. 2022. Transcriptomic Analyses of Grapevine Leafroll-Associated Virus 3 Infection in Leaves and Berries of 'Cabernet Franc'. *Viruses* 14(8):1831. | 10.3390/v14081831 / 36016453 | RNA-Seq; Cabernet Franc leaves + berries, 3 stages | 6 contrasts. Pre-filtered to FDR < 0.05. |
| **Powdery mildew** (biotroph fungus) | Fung RWM, Gonzalo M, Fekete C, Kovacs LG, He Y, Marsh E, McIntyre LM, Schachtman DP, Qiu W. 2008. Powdery mildew induces defense-oriented reprogramming of the transcriptome in a susceptible but not in a resistant grapevine. *Plant Physiol* 146(1):236–249. | 10.1104/pp.107.108712 / 17993546 | Affymetrix Vitis GeneChip; Cabernet Sauvignon leaves | ID space: Affymetrix probeset + Arabidopsis homolog (AT number). 0.5% coverage against our crosswalk (known limitation — Tier 2 bridge needed). |
| **Downy mildew** (oomycete biotroph) | Polesani M, Bortesi L, Ferrarini A, Zamboni A, Fasoli M, Zadra C, Lovato A, Pezzotti M, Delledonne M, Polverari A. 2010. General and species-specific transcriptional responses to downy mildew infection in a susceptible (*Vitis vinifera*) and a resistant (*V. riparia*) grapevine species. *BMC Genomics* 11:117. | 10.1186/1471-2164-11-117 / 20167053 | Combimatrix microarray; V. vinifera leaf disks 24 hpi | ID space: TIGR Tentative Consensus (TC#). 0% coverage — sequence-based TC→VIT mapping needed for Tier 2. |
| **Botrytis cinerea** (necrotroph on ripening berries) | Agudelo-Romero P, Erban A, Rego C, Carbonell-Bejerano P, Nascimento T, Sousa L, Martínez-Zapater JM, Kopka J, Fortes AM. 2015. Transcriptome and metabolome reprogramming in *Vitis vinifera* cv. Trincadeira berries upon infection with *Botrytis cinerea*. *J Exp Bot* 66(7):1769–1785. | 10.1093/jxb/eru517 / 25675955 | NimbleGen microarray; Trincadeira berries, EL33 + EL35 | 2 contrasts. `12X V1 ID` = VIT_*. 60% coverage. |
| **Heat stress** (abiotic) | Liu GT, Wang JF, Cramer G, Dai ZW, Duan W, Xu HG, Wu BH, Fan PG, Wang LJ, Li SH. 2012. Transcriptomic analysis of grape (*Vitis vinifera* L.) leaves during and after recovery from heat stress. *BMC Plant Biol* 12:174. | 10.1186/1471-2229-12-174 / 23016701 | Affymetrix Grape Genome microarray; Cabernet Sauvignon leaves | Same chip as Fung 2008. ID space: probeset + GenBank EST accession. 0% coverage — shared probe→VIT table with Fung unlocks both. |

Coverage summary produced by rule 11 (see `results/11_harmonized/coverage_report.csv`):

| Paper | Rows mapped | Coverage |
|---|---|---|
| Song 2022 | 304–2088 per contrast | 68–91 % |
| Blanco-Ulate 2017 | 17,846 per contrast | 68.7 % |
| Agudelo-Romero 2015 | 5,516 per contrast | 60.2 % |
| Fung 2008 | 2 | 0.5 % |
| Liu 2012 | 0 | 0 % |
| Polesani 2010 | 0 | 0 % |

The three microarray-era papers with probeset / TC identifiers contribute to
the comparative analysis only after a Tier 2 ID-space bridge is built (see
Section 4).

---

## 2. Supplementary files ingested

| File in `data/literature/` | Source paper | DEG content |
|---|---|---|
| `erw506_suppl_Supplementary_Table_S5.xlsx` | Blanco-Ulate 2017 | Table S5, full transcriptome, 25,997 VIT genes × 8 contrasts |
| `Song et al. RNA-Seq grapevine-GLRaV-3 STable S4 ...xlsx` | Song 2022 | Table S4, 6 sheets (tissue × stage), pre-filtered FDR<0.05 |
| `108712Table_S1.xls` | Fung 2008 | Table S1, 626 PM-responsive probes, time-course |
| `12864_2009_2711_MOESM6_ESM.XLS` | Polesani 2010 | Additional File 6, 4 side-by-side subtables (All / Vv-specific / Vr-specific / common) |
| `jexbot132548_file001.xlsx` | Agudelo-Romero 2015 | Sheet `Table S4SignificantGenes2.0fold`, 9,168 probes × 4 contrasts |
| `12870_2012_1124_MOESM7_ESM.docx` | Liu 2012 | Additional File 7, 490-row Word table (heat + recovery FC) |

Ingestion is one-shot: `python3 scripts/ingest_tier1.py` reads the six
source files, emits 18 tidy CSVs + 18 YAML sidecars into
`data/literature/`, one per contrast, following the naming convention
below.

---

## 3. File schema

### Naming convention

`<firstauthor><year>_<stress>_<cultivar>_<tissue>_<stage>.csv`

Examples:
- `BlancoUlate2017_GRBV_Zinfandel_berry_healdsburg_harvest.csv`
- `Song2022_GLRaV3_CabernetFranc_leaf_el38.csv`
- `AgudeloRomero2015_Botrytis_Trincadeira_berry_el35.csv`
- `Fung2008_PowderyMildew_CabernetSauvignon_leaf_24h.csv`

### Required CSV columns

| Column | Type | Description |
|---|---|---|
| `gene_id` | str | identifier as published (VIT_\*, Vitvi\*, LOC\*, Symbol, or AT#) |
| `log2FC` | float | log₂ fold change, stress vs control (positive = up in stress) |
| `pvalue` | float | raw p-value from the DE test (NaN if paper only reports padj) |
| `padj` | float | FDR-adjusted p-value (set to threshold if paper was pre-filtered) |

### Optional columns

| Column | Default | Description |
|---|---|---|
| `id_type` | inferred | `VIT`, `Vvi`, `LOC`, or `Symbol` |
| `tested` | 1 | `1` if in the test universe, `0` otherwise |
| `baseMean` | NaN | mean expression (informational) |

### YAML sidecar

Each CSV has a `<basename>.yaml` with metadata:

```yaml
stress: "GRBV"             # GRBV | GLRaV-3 | PowderyMildew | DownyMildew | Botrytis | Heat | ...
cultivar: "Zinfandel (Healdsburg)"
tissue: "berry"            # leaf | berry | rootstock | xylem | phloem
stage: "Harvest"
reference: "Blanco-Ulate et al. 2017, J Exp Bot 68(5):1225-1238, doi:10.1093/jxb/erw506"
platform: "RNA-Seq"
deg_threshold_padj: 0.05
deg_threshold_lfc: 1.0
full_panel: true           # true if the CSV contains all tested genes; false if pre-filtered
note: "Optional free-text caveat."
```

The pipeline uses `stress` to group datasets for the GRBV-unique /
viral-shared / pathogen-generic / stress-generic classifier. Misspelling it
will silently shunt the dataset into the "Mixed" bucket.

---

## 4. Two-tier analysis plan

### Tier 1 — List-level meta-analysis (runs today)

1. Extract published DEG tables from each paper's supplementary files
   (done, six files in this directory).
2. Harmonize gene IDs to NCBI LOC\* via the ID crosswalk at
   `resources/id_crosswalk.tsv`. The crosswalk (24,817 rows) combines:
   - NCBI current-assembly records (Vvi ↔ LOC ↔ Symbol) from
     `All_Plants.gene_info.gz`;
   - Ensembl Plants release-50 12X.v1 records (VIT\_\*s\*g\* ↔ LOC) from
     `Vitis_vinifera.12X.50.entrez.tsv.gz`.
3. Compute for every pair of signatures:
   - Fisher exact test on overlap (odds ratio, p, BH FDR)
   - Spearman correlation of log₂FC on the intersection of tested genes
   - Signed concordance (fraction of overlapping DEGs with same direction)
4. Classify each GRBV DEG into one of four categories:
   - **GRBV-unique**: significant in GRBV only
   - **GRBV + viral-shared**: GRBV + ≥ 1 other virus, not in mildews / heat
   - **Generic pathogen**: GRBV + ≥ 1 mildew
   - **Generic stress**: significant in ≥ 3 datasets including non-pathogen

Tier 1 outputs (current run with 6 papers, 55 GRBV DEGs at padj<0.05 +
|logFC|≥1):

| Output | Location |
|---|---|
| Harmonized per-dataset tables | `results/11_harmonized/*.harmonized.csv` |
| Coverage report | `results/11_harmonized/coverage_report.csv` |
| Per-pair overlap statistics | `results/12_metaanalysis/overlap_matrix.csv` |
| log₂FC correlation matrix | `results/12_metaanalysis/fc_correlation_matrix.csv` |
| Per-gene specificity classes | `results/12_metaanalysis/specificity_classes.csv` |
| Specificity summary | `results/12_metaanalysis/specificity_summary.csv` |
| UpSet plot, FC heatmap, bar plot | `results/12_metaanalysis/*.png` |

Headline results (2026-04-23 run):

- 55 GRBV DEGs in our study; 35 GRBV-unique, 20 viral-shared.
- **Strongest cross-study replication: GRBV (Merlot, this study) vs
  Blanco-Ulate 2017 (Zinfandel, Oakville, Veraison)**: 4 of 39 shared DEGs,
  Fisher p = 2.0 × 10⁻⁸, odds ratio = 225, Spearman ρ(log₂FC) = 0.51. Same
  pair at Oakville Harvest: 4/39, ρ = 0.47.
- GLRaV-3 × Cabernet Franc (Song 2022) shows weak or no overlap with our
  GRBV signature (ρ ≈ 0 at most stages), consistent with the two viruses
  producing distinct host-transcriptome signatures.
- Agudelo-Romero Botrytis contributes to the overlap network but does not
  drive the GRBV-unique / viral-shared split at the current thresholds.

### Tier 2 — Effect-size / raw-count re-analysis (future)

For the three microarray-era papers with unresolved ID spaces (Fung 2008,
Polesani 2010, Liu 2012), Tier 2 requires one or more of the following
bridges before their DEG signatures contribute meaningfully:

1. **Affymetrix Vitis GeneChip probeset → VIT** (unlocks Fung 2008 + Liu
   2012; same chip). Two feasible routes: (a) download the Vitis
   GeneChip annotation from NetAffx / Bioconductor `vitis.db` if still
   available; (b) BLAST probe sequences against V. vinifera CDS.
2. **TIGR Tentative Consensus (TC#) → VIT** (unlocks Polesani 2010).
   TIGR's Grapevine Gene Index is deprecated; the mapping typically
   comes from the Polesani authors' design file or from re-BLASTing TC
   consensus sequences.
3. **Raw-count re-analysis from GEO / SRA** for all six papers using a
   single DESeq2 pipeline. Most robust option, 1–2 weeks of work.

Until the bridges exist, the current Tier 1 output treats these three
studies as "tested but zero overlap" in the overlap matrix (conservative
but correct).

---

## 5. Rebuild procedure

From the repo root:

```bash
# Ingest the six supplementary files -> 18 tidy CSVs + YAML sidecars
python3 scripts/ingest_tier1.py

# Harmonize to NCBI LOC_ IDs
cd farm_pipeline
python3 scripts/11_harmonize_literature_degs.py \
    --literature_dir data/literature \
    --crosswalk      resources/id_crosswalk.tsv \
    --out_dir        results/11_harmonized \
    --out_coverage   results/11_harmonized/coverage_report.csv

# Run the comparative meta-analysis
python3 scripts/12_comparative_metaanalysis.py \
    --counts         ~/Data_Analysis/counts_matrix.txt \
    --meta           ~/Data_Analysis/sample_metadata.txt \
    --harmonized_dir results/11_harmonized \
    --coverage       results/11_harmonized/coverage_report.csv \
    --out_dir        results/12_metaanalysis \
    --padj_thresh    0.05 \
    --lfc_thresh     1.0
```

To rebuild the crosswalk (only needed if new papers are added with
different ID formats):

```bash
# Current NCBI annotation side (Vvi, LOC, Symbol)
curl -sSL https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Plants/All_Plants.gene_info.gz \
  -o /tmp/all_plants.gene_info.gz
# Ensembl 12X side (VIT_*s*g* -> LOC)
curl -sSL https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-50/tsv/vitis_vinifera/Vitis_vinifera.12X.50.entrez.tsv.gz \
  -o /tmp/Vv_12x_entrez.tsv.gz
# Merge into resources/id_crosswalk.tsv (see scripts/build_id_crosswalk.sh).
```
