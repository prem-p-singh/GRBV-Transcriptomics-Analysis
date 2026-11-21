#!/bin/bash
#SBATCH --job-name=grbv_perm
#SBATCH --account=publicgrp          # free tier
#SBATCH --partition=low              # preemptible; requeued tasks skip finished draws
#SBATCH --time=01:00:00              # per task
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16           # dream workers for one draw
#SBATCH --mem=48G
#SBATCH --requeue
#SBATCH --output=logs/perm-%A_%a.out
#SBATCH --error=logs/perm-%A_%a.err

source /etc/profile.d/modules.sh
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
env/bin/Rscript 07d_perm_task.R
