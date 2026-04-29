#!/bin/bash
# =============================================================================
# run_pipeline_cactus.sh
# Submit the cactus_em pipeline as a chained Slurm dependency series.
#   00_setup_cactus_clones ── 02_run_hack_cactus ── 03_run_shorts_cactus ── 05c_cactus_em
# =============================================================================
set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE=${1:-"${SCRIPTS_DIR}/config_cactus_em.sh"}
mkdir -p "${SCRIPTS_DIR}/../logs"

JOB0=$(sbatch --parsable "${SCRIPTS_DIR}/00_setup_cactus_clones.sh" "${CONFIG_FILE}")
echo "Submitted 00_setup_cactus_clones.sh   as job ${JOB0}"

JOB2=$(sbatch --parsable --dependency=afterok:${JOB0} \
       "${SCRIPTS_DIR}/02_run_hack_cactus.sh" "${CONFIG_FILE}")
echo "Submitted 02_run_hack_cactus.sh        as job ${JOB2} (afterok:${JOB0})"

JOB3=$(sbatch --parsable --dependency=afterok:${JOB2} \
       "${SCRIPTS_DIR}/03_run_shorts_cactus.sh" "${CONFIG_FILE}")
echo "Submitted 03_run_shorts_cactus.sh      as job ${JOB3} (afterok:${JOB2})"

JOB5=$(sbatch --parsable --dependency=afterok:${JOB3} \
       "${SCRIPTS_DIR}/05c_cactus_em.sh" "${CONFIG_FILE}")
echo "Submitted 05c_cactus_em.sh             as job ${JOB5} (afterok:${JOB3})"

echo ""
echo "Job chain: ${JOB0} → ${JOB2} → ${JOB3} → ${JOB5}"
echo "Run eval after JOB5 completes:"
echo "  python ${SCRIPTS_DIR}/06_eval_cactus_em.py --config ${CONFIG_FILE} --size 1kb"
