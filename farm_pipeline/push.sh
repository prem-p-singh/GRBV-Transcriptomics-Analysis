#!/bin/bash
# Push pipeline + inputs from Mac to FARM.
# Run from ~/Data_Analysis/farm_pipeline/ :  bash push.sh
set -euo pipefail

REMOTE="farm:~/grbv_extended/"

# Stage raw inputs into data/ (if not already there)
DATA_DIR="data"
SRC="$HOME/Data_Analysis"
mkdir -p "$DATA_DIR"
cp -f "$SRC/counts_matrix.txt"     "$DATA_DIR/"
cp -f "$SRC/sample_metadata.txt"   "$DATA_DIR/"
cp -f "$SRC/Virus titer.txt"       "$DATA_DIR/Virus_titer.txt"   # rename (no space)
cp -f "$SRC/MEs.csv"               "$DATA_DIR/"
cp -f "$SRC/Module_genes.xlsx"     "$DATA_DIR/"
cp -f "$SRC/daily_2021.csv"        "$DATA_DIR/"
cp -f "$SRC/daily_2022.csv"        "$DATA_DIR/"

# Rsync the whole pipeline tree (exclude results/logs so we don't overwrite remote output)
rsync -avz --progress \
    --exclude='results/' \
    --exclude='logs/' \
    --exclude='resources/genome/' \
    --exclude='resources/motifs/' \
    --exclude='resources/*.txt' \
    --exclude='__pycache__/' \
    --exclude='.DS_Store' \
    ./ "$REMOTE"

echo "Pushed to $REMOTE"
