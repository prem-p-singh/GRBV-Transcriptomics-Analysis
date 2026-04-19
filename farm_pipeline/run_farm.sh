#!/bin/bash
# Launch the Snakemake pipeline on FARM. Run from ~/grbv_extended/ on FARM.
set -euo pipefail

set +u
source /etc/profile.d/modules.sh
module load conda/base
conda activate grbv-extended
set -u

snakemake --profile profiles/farm/ --configfile config/config.yaml "$@"
