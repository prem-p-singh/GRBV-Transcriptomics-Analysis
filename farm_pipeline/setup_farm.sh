#!/bin/bash
# FARM bootstrap: Python env (conda) + R packages (LMOD R/4.4.2 + user lib) + resources.
# Run from ~/grbv_extended/ on FARM:  bash setup_farm.sh
set -euo pipefail

set +u
source /etc/profile.d/modules.sh
module load conda/base
set -u

ENV_NAME="grbv-extended"
RES="resources"
RLIBS="$HOME/grbv_extended/Rlibs"

mkdir -p "$RES/genome" "$RES/motifs" "$RLIBS" logs

# ---------- 1. Python conda env (lean, no R) ----------
if [ ! -x "$HOME/.conda/envs/$ENV_NAME/bin/python" ]; then
    echo "[setup] creating conda env $ENV_NAME (python + tooling only)"
    [ -d "$HOME/.conda/envs/$ENV_NAME" ] && rm -rf "$HOME/.conda/envs/$ENV_NAME"
    conda env create -f envs/grbv-extended.yaml
else
    echo "[setup] conda env $ENV_NAME already exists"
fi

# ---------- 2. R via LMOD + user library ----------
echo "[setup] loading R/4.4.2 LMOD module"
set +u; module load R/4.4.2; set -u
export GRBV_RLIBS="$RLIBS"
echo "[setup] installing R packages into $RLIBS (this takes a while)"
Rscript envs/install_r_packages.R

# ---------- 3. Reference genome + GFF ----------
conda activate "$ENV_NAME"
ACC=$(awk '/ncbi_accession/{print $2}' config/config.yaml)
if [ ! -f "$RES/genome/genome.fna" ]; then
    echo "[setup] downloading NCBI genome $ACC"
    cd "$RES/genome"
    datasets download genome accession "$ACC" --include genome,gff3
    unzip -o ncbi_dataset.zip
    cp ncbi_dataset/data/"$ACC"/*_genomic.fna genome.fna
    cp ncbi_dataset/data/"$ACC"/genomic.gff    genome.gff
    rm -rf ncbi_dataset.zip README.md md5sum.txt ncbi_dataset
    cd - >/dev/null
else
    echo "[setup] genome already present"
fi

# ---------- 4. PlantTFDB Vvi TF list ----------
if [ ! -f "$RES/Vvi_TF_list.txt" ]; then
    echo "[setup] downloading PlantTFDB Vvi TF list"
    curl -sSL -o "$RES/Vvi_TF_list.txt.gz" \
        "http://planttfdb.gao-lab.org/download/TF_list/Vvi_TF_list.txt.gz"
    gunzip -f "$RES/Vvi_TF_list.txt.gz"
else
    echo "[setup] TF list already present"
fi

# ---------- 5. JASPAR plants motif DB ----------
if [ ! -f "$RES/motifs/jaspar_plants.meme" ]; then
    echo "[setup] downloading JASPAR2024 CORE plants motifs"
    curl -sSL -o "$RES/motifs/jaspar_plants.meme" \
        "https://jaspar.genereg.net/download/data/2024/CORE/JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt"
else
    echo "[setup] motif DB already present"
fi

# ---------- 6. ID crosswalk for comparative meta-analysis (rules 11-12) ----------
if [ ! -f "$RES/id_crosswalk.tsv" ]; then
    echo "[setup] building VIT/Vvi/Symbol -> LOC crosswalk from NCBI GFF"
    python scripts/_build_id_crosswalk.py \
        --gff "$RES/genome/genome.gff" \
        ${ORTHOLOG_MAP:+--ortholog_map "$ORTHOLOG_MAP"} \
        --out "$RES/id_crosswalk.tsv"
else
    echo "[setup] id_crosswalk.tsv already present"
fi

# ---------- 7. Stub literature directory if empty ----------
mkdir -p data/literature
if [ -z "$(ls -A data/literature/*.csv 2>/dev/null)" ]; then
    echo "[setup] data/literature/ is empty -- rule 11 will skip until you populate it."
    echo "        See docs/LITERATURE_SCHEMA.md for the required CSV format."
fi

echo "[setup] disk usage:"
du -sh "$RES"/* 2>/dev/null || true
du -sh "$RLIBS" 2>/dev/null || true
df -h ~ | tail -1
echo "[setup] done"
