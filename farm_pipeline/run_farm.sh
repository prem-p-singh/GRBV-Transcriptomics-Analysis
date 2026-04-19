#!/bin/bash
# Launch the Snakemake pipeline on FARM. Run from ~/grbv_extended/ on FARM.
set -euo pipefail

set +u
source /etc/profile.d/modules.sh
module load conda/base
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate grbv-extended
module load R/4.4.2
export GRBV_RLIBS="$HOME/grbv_extended/Rlibs"
set -u

snakemake --profile profiles/farm/ --configfile config/config.yaml "$@"
