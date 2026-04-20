# `data/literature/` — published Vitis DEG tables for the comparative meta-analysis

Rules 11 and 12 of the FARM pipeline consume published DEG tables from this
directory. Each file is one comparison (one stress vs its own control) from
one paper. The pipeline maps heterogeneous identifiers to NCBI `LOC_` IDs,
computes overlap statistics against the GRBV DEG signature, and classifies
each GRBV DEG as GRBV-unique / viral-shared / pathogen-generic / stress-
generic.

## File naming convention

`<firstauthor><year>_<stress>_<cultivar>_<tissue>.csv`

Examples:
- `BlancoUlate2017_GLRaV3_Cabernet_leaf.csv`
- `Bertazzon2020_GPGV_Pinot_leaf.csv`
- `Amrine2015_PM_Cabernet_leaf.csv`
- `Perazzolli2012_DM_Pinot_leaf.csv`
- `Haider2017_Heat_Thompson_leaf.csv`

## Required CSV schema

Every file must have the following columns (other columns are ignored):

| Column | Type | Description |
|---|---|---|
| `gene_id` | str | identifier as published (VIT_\*, Vitvi\*, LOC\*, or symbol) |
| `log2FC` | float | log₂ fold change, stress vs control, sign: up in stress = positive |
| `pvalue` | float | raw p-value from the DE test |
| `padj` | float | FDR-adjusted p-value |

Optional columns the pipeline respects:

| Column | Type | Description |
|---|---|---|
| `id_type` | str | `VIT`, `Vvi`, `LOC`, or `Symbol`. Inferred from the prefix if missing. |
| `tested` | int | `1` if the gene was in the test universe, `0` otherwise (default: 1) |
| `baseMean` | float | mean expression — purely informational |

## Optional sidecar YAML

Place a `<basename>.yaml` next to the CSV to attach metadata:

```yaml
stress: "GLRaV-3"          # GLRaV-3 | GPGV | GRBV | PowderyMildew | DownyMildew | Heat | Drought | Cold | ...
cultivar: "Cabernet Sauvignon"
tissue: "leaf"             # leaf | berry | rootstock | xylem | phloem
reference: "Blanco-Ulate et al. 2017, Plant Physiol 174:1054-1068, doi:10.1104/pp.17.00420"
deg_threshold_padj: 0.05
deg_threshold_lfc: 1
notes: "Symptomatic late-stage samples only; healthy control from same vineyard."
```

The pipeline uses the `stress` field to group datasets (viral vs mildew vs abiotic) when classifying GRBV DEG specificity. Misspelling the stress group will silently put the dataset into "Mixed".

## Candidate datasets to ingest

Five datasets give good coverage of the Vitis stress landscape:

| Stress | Suggested study | Accession |
|---|---|---|
| GLRaV-3 | Blanco-Ulate et al. 2017 *Plant Physiol* | GEO GSE98474 |
| GPGV | Bertazzon et al. 2020 *Viruses* | SRA PRJNA574257 |
| Powdery mildew (*E. necator*) | Amrine et al. 2015 *BMC Plant Biol* | GEO GSE63859 |
| Downy mildew (*P. viticola*) | Perazzolli et al. 2012 *BMC Genomics* | GEO GSE34784 |
| Heat stress (abiotic control) | Haider et al. 2017 *BMC Plant Biol* | GEO GSE114322 |

For each study, download the supplementary DEG table (Excel or CSV), reduce
to the four required columns, save to this directory as a CSV, and write
the YAML sidecar. See the template below.

## Minimal template

```csv
gene_id,log2FC,pvalue,padj
VIT_12s0134g00070,1.42,3.2e-05,7.8e-04
VIT_04s0023g01030,-2.05,1.1e-08,3.4e-06
LOC100247440,0.91,8.3e-04,1.2e-02
```

```yaml
# dataset.yaml
stress: "GLRaV-3"
cultivar: "Cabernet Sauvignon"
tissue: "leaf"
reference: "Example et al. 2020"
deg_threshold_padj: 0.05
deg_threshold_lfc: 1
```

## ID crosswalk

Rule 11 uses `resources/id_crosswalk.tsv` to map non-`LOC_` identifiers to
`LOC_` IDs. It is built once during setup by `scripts/_build_id_crosswalk.py`
from the NCBI GFF (which carries `old_locus_tag=VIT_*` for most genes).

For higher coverage on Vitvi* (12X.v2 / PN40024 annotation) IDs, download
an ortholog map from [Ensembl Plants BioMart](http://plants.ensembl.org/biomart)
with columns `Gene stable ID` (Vitvi), `VIT_gene_stable_ID`, `NCBI gene (LOC) ID`,
`Gene name`, rename to `Vvi_id`, `VIT_id`, `LOC_id`, `Symbol`, and pass the
file path via the `--ortholog_map` flag of the helper script. Without this
extra map, studies using Vvi identifiers will suffer ~50 % mapping loss.

## Expected outputs once rule 11 runs

- `results/11_harmonized/<dataset>.harmonized.csv` — per dataset, `LOC_id` + DE stats
- `results/11_harmonized/coverage_report.csv` — per-dataset mapping success

If `coverage_pct` < 70 %, the dataset is probably using an annotation the
crosswalk doesn't cover; add the appropriate ortholog table and rerun.
