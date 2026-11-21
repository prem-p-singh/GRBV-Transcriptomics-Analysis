#!/bin/bash
#SBATCH --job-name=grbv_env
#SBATCH --account=publicgrp          # free tier; ppssingh's only account
#SBATCH --partition=low              # preemptible; a rebuild just restarts
#SBATCH --time=02:00:00              # wall-clock limit
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --requeue
#SBATCH --output=logs/env-%j.out
#SBATCH --error=logs/env-%j.err

# Builds an R environment. Bioconda has no 1.38.1 (the local version), so 1.40.2 on
# R 4.5 is used; 07d checks that it reproduces the primary counts before the array.
# Package caches stay in per-job /tmp, so only the finished env lands in home.
source /etc/profile.d/modules.sh
module load conda
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export CONDA_PKGS_DIRS="$TMPDIR/conda_pkgs"
rm -rf env
conda create -y -p "$PWD/env" --override-channels -c conda-forge -c bioconda \
  "bioconductor-variancepartition=1.40.2" r-data.table
env/bin/Rscript -e 'for (p in c("variancePartition","limma","lme4","data.table")) cat(p, as.character(packageVersion(p)), "\n"); cat(R.version.string, "\n")'
du -sh env
