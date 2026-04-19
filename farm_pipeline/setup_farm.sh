#!/bin/bash
# One-shot FARM bootstrap: create conda env + download resources.
# Run from ~/grbv_extended/ on FARM:  bash setup_farm.sh
set -euo pipefail

set +u
source /etc/profile.d/modules.sh
module load conda/base
set -u

ENV_NAME="grbv-extended"
RES="resources"

# ---------- 1. conda env ----------
if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "[setup] creating conda env $ENV_NAME"
    conda env create -f envs/grbv-extended.yaml
else
    echo "[setup] conda env $ENV_NAME already exists"
fi
conda activate "$ENV_NAME"

mkdir -p "$RES/genome" "$RES/motifs" logs

# ---------- 2. Vitis vinifera reference genome + GFF ----------
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

# ---------- 3. PlantTFDB TF list for V. vinifera ----------
if [ ! -f "$RES/Vvi_TF_list.txt" ]; then
    echo "[setup] downloading PlantTFDB Vvi TF list"
    curl -sSL -o "$RES/Vvi_TF_list.txt.gz" \
        "http://planttfdb.gao-lab.org/download/TF_list/Vvi_TF_list.txt.gz"
    gunzip -f "$RES/Vvi_TF_list.txt.gz"
else
    echo "[setup] TF list already present"
fi

# ---------- 4. JASPAR plants motif database (MEME format) ----------
if [ ! -f "$RES/motifs/jaspar_plants.meme" ]; then
    echo "[setup] downloading JASPAR2024 CORE plants motifs"
    curl -sSL -o "$RES/motifs/jaspar_plants.meme" \
        "https://jaspar.genereg.net/download/data/2024/CORE/JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt"
else
    echo "[setup] motif DB already present"
fi

echo "[setup] done. disk usage:"
du -sh "$RES"/* 2>/dev/null
df -h ~ | tail -1
