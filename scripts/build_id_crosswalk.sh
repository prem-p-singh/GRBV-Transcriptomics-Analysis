#!/usr/bin/env bash
# Build farm_pipeline/resources/id_crosswalk.tsv: a VIT_*s*g* (12X.v1) and
# Vvi (current Ensembl Plants) -> NCBI LOC_ crosswalk. Merges two sources:
#
#   1. NCBI All_Plants.gene_info.gz  -> current Vvi + LOC + Symbol
#   2. Ensembl Plants release-50 entrez TSV -> legacy VIT_*s*g* -> LOC
#
# Both public, both free. Output is ~25k rows. Re-run only when adding new
# literature papers with identifiers not already covered.
set -euo pipefail

RES="$(dirname "$0")/../farm_pipeline/resources"
mkdir -p "$RES"

echo "[1/4] downloading NCBI All_Plants.gene_info.gz"
curl -sSL https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Plants/All_Plants.gene_info.gz \
    -o /tmp/all_plants.gene_info.gz

echo "[2/4] downloading Ensembl Plants 12X entrez TSV"
curl -sSL https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-50/tsv/vitis_vinifera/Vitis_vinifera.12X.50.entrez.tsv.gz \
    -o /tmp/Vv_12x_entrez.tsv.gz

echo "[3/4] extracting V. vinifera current-annotation side"
gunzip -c /tmp/all_plants.gene_info.gz | awk -F'\t' -v OFS='\t' '
BEGIN { print "VIT_id","Vvi_id","LOC_id","Symbol","source" }
$1=="29760" {
    # Skip VIT col (new-format NCBI tags conflict with 12X.v1); keep Vvi from dbXrefs
    vvi = ""
    if (match($6, /Vitis[0-9]+g[0-9]+/)) vvi = substr($6, RSTART, RLENGTH)
    sym = ($5 != "-" && $5 != "") ? $5 : ""
    if (vvi != "" || sym != "") print "", vvi, $3, sym, "ncbi_gene_info"
}' > /tmp/cw_current.tsv

echo "[4/4] extracting Ensembl 12X VIT_*s*g* -> LOC, merging"
( echo -e "VIT_id\tVvi_id\tLOC_id\tSymbol\tsource"
  gunzip -c /tmp/Vv_12x_entrez.tsv.gz | awk -F'\t' -v OFS='\t' '
    NR==1 { next }
    $1 ~ /^VIT_/ && $4 ~ /^[0-9]+$/ { print $1, "", "LOC"$4, "", "ensembl_12x_entrez" }
  ' | sort -u ) > /tmp/cw_12x.tsv

python3 - <<'PY'
import pandas as pd
from pathlib import Path

OUT = Path("/Users/prempratapsingh/Data_Analysis/farm_pipeline/resources/id_crosswalk.tsv")
a = pd.read_csv("/tmp/cw_current.tsv", sep="\t", dtype=str).fillna("")
b = pd.read_csv("/tmp/cw_12x.tsv",     sep="\t", dtype=str).fillna("")

# Enrich 12X rows with Vvi + Symbol from current NCBI by joining on LOC
a_by_loc = a[a["LOC_id"].ne("")].drop_duplicates("LOC_id", keep="first")[["LOC_id","Vvi_id","Symbol"]]
b_enriched = b.merge(a_by_loc, on="LOC_id", how="left", suffixes=("","_cur")).fillna("")
b_enriched["Vvi_id"] = b_enriched["Vvi_id_cur"]
b_enriched["Symbol"] = b_enriched["Symbol_cur"]
b_enriched = b_enriched[["VIT_id","Vvi_id","LOC_id","Symbol","source"]]

# Append NCBI rows whose LOC is not already covered by 12X mapping
locs = set(b_enriched["LOC_id"])
a_extra = a[~a["LOC_id"].isin(locs)].copy()
a_extra["VIT_id"] = ""
merged = pd.concat([b_enriched, a_extra], ignore_index=True).drop_duplicates(
    subset=["VIT_id","Vvi_id","LOC_id","Symbol"], keep="first")
merged.to_csv(OUT, sep="\t", index=False)
print(f"wrote {OUT} ({len(merged)} rows)")
PY

echo "done."
