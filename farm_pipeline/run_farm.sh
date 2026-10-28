#!/bin/bash
# Launch the Snakemake pipeline on FARM. Run from ~/grbv_extended/ on FARM.
set -euo pipefail

set +u
source /etc/profile.d/modules.sh

# The Snakemake controller runs as its own job, not on the login node. Start this inside
# tmux: srun stops if the ssh session drops. `low` is preemptible; rerun to resume
# (rerun-incomplete is on in the profile).
if [ -z "${SLURM_JOB_ID:-}" ]; then
    exec srun -A publicgrp -p low -t "${RUN_TIME:-3-00:00:00}" -c "${RUN_CPUS:-1}" \
        --mem="${RUN_MEM:-4G}" bash "$0" "$@"
fi
# This job's SLURM_*/SBATCH_* variables would otherwise be exported into every job
# Snakemake submits (e.g. SLURM_MEM_PER_CPU clashing with the rule's --mem).
# SLURM_CONF is kept so sbatch still finds the cluster config.
for v in $(compgen -e | grep -E '^(SLURM_|SBATCH_|SRUN_|SALLOC_)' | grep -vx SLURM_CONF); do
    unset "$v"
done

module load conda/base
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate grbv-extended
module load R/4.4.2
export GRBV_RLIBS="$HOME/grbv_extended/Rlibs"
set -u

snakemake --profile profiles/farm/ --configfile config/config.yaml "$@"
