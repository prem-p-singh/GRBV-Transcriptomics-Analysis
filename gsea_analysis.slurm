#!/bin/bash
#SBATCH -D /home/ppssingh/gsea_run
#SBATCH -o /home/ppssingh/gsea_run/gsea_stdout-%j.txt
#SBATCH -e /home/ppssingh/gsea_run/gsea_stderr-%j.txt
#SBATCH -J gsea_grbv
#SBATCH -p low
#SBATCH -A publicgrp
#SBATCH -t 4:00:00
#SBATCH --mem 16G
#SBATCH -n 4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ppssingh@ucdavis.edu

start=`date +%s`
echo $HOSTNAME

outpath="gsea_output"
mkdir -p ${outpath}

module load python/3.11

echo "Python version:"
python3 --version

echo "Checking dependencies..."
python3 -c "import numpy, pandas, scipy, matplotlib, seaborn, goatools; print('All dependencies OK')"

echo "============================================"
echo "Running GSEA (vectorized, 500 permutations)"
echo "============================================"

call="python3 gsea_analysis.py"

echo $call
eval $call

echo "============================================"
echo "Checking outputs..."
echo "============================================"

for f in gsea_results_all_conditions.csv 37_gsea_enrichment_plots.png 38_gsea_vs_fisher_concordance.png 39_gsea_nes_heatmap.png; do
    if [[ -f $f ]]; then
        echo "  OK: $f ($(ls -lh $f | awk '{print $5}'))"
        cp $f ${outpath}/
    else
        echo "  MISSING: $f"
    fi
done

end=`date +%s`
runtime=$((end-start))
echo $runtime
