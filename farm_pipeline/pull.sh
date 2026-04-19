#!/bin/bash
# Pull results + logs from FARM back to Mac.
set -euo pipefail

rsync -avz --progress \
    farm:~/grbv_extended/results/ ./results/

rsync -avz --progress \
    farm:~/grbv_extended/logs/ ./logs/

echo "Pulled results/ and logs/"
