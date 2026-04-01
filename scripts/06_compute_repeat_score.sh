#!/bin/bash
#SBATCH --job-name=repeat_score
#SBATCH --output=logs/repeat_score_%j.out
#SBATCH --error=logs/repeat_score_%j.err
#SBATCH --time=00:15:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# =============================================================================
# 06_compute_repeat_score.sh
# Compute genome-wide k-mer repeat score for a registered position and
# update data/positions_registry.tsv in-place.
#
# Runs in parallel with the main pipeline (no dependencies).
# Only updates the registry for rep-style labels (rep1, rep2, …).
# Legacy pos-style labels (pos10mb, pos20mb, …) are skipped.
# =============================================================================
set -euo pipefail

CONFIG_FILE=${1:?"Usage: $0 <config_file>"}
source "${CONFIG_FILE}"

# Only act on rep-style labels
if [[ "${POS_LABEL}" != rep* ]]; then
    echo "POS_LABEL=${POS_LABEL} is not a rep label — skipping repeat score computation."
    exit 0
fi

echo "Computing repeat score for ${POS_LABEL} at pos ${SV_START_0} …"

python3 "${WORK}/scripts/compute_repeat_score.py" \
    --rep-id "${POS_LABEL}" \
    --pos-bp "${SV_START_0}"
