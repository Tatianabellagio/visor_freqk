#!/bin/bash
# =============================================================================
# run_pipeline.sh
# Purpose: Submit the full VISOR + freqk pipeline (steps 01–05) as a chained
#          series of Slurm jobs for a given config file (DEL or INS).
#
# Usage:
#   cd /home/tbellagio/scratch/pang/visor_freqk
#   bash scripts/run_pipeline.sh scripts/config_sv_deletions.sh
#   bash scripts/run_pipeline.sh scripts/config_sv_insertions.sh
# =============================================================================
set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE=${1:-"${SCRIPTS_DIR}/config_sv_deletions.sh"}

if [[ ! -f "${CONFIG_FILE}" ]]; then
  echo "Config file not found: ${CONFIG_FILE}" >&2
  exit 1
fi

echo "Using config: ${CONFIG_FILE}"

# Submit step 01
JOB1=$(sbatch --parsable "${SCRIPTS_DIR}/01_make_beds.sh" "${CONFIG_FILE}")
echo "Submitted 01_make_beds.sh as job ${JOB1}"

# Step 02 depends on 01
JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} "${SCRIPTS_DIR}/02_run_hack.sh" "${CONFIG_FILE}")
echo "Submitted 02_run_hack.sh as job ${JOB2} (afterok:${JOB1})"

# Step 03 depends on 02
JOB3=$(sbatch --parsable --dependency=afterok:${JOB2} "${SCRIPTS_DIR}/03_run_shorts.sh" "${CONFIG_FILE}")
echo "Submitted 03_run_shorts.sh as job ${JOB3} (afterok:${JOB2})"

# Step 04 depends on 03
JOB4=$(sbatch --parsable --dependency=afterok:${JOB3} "${SCRIPTS_DIR}/04_make_vcf.sh" "${CONFIG_FILE}")
echo "Submitted 04_make_vcf.sh as job ${JOB4} (afterok:${JOB3})"

# Step 05 depends on 04
JOB5=$(sbatch --parsable --dependency=afterok:${JOB4} "${SCRIPTS_DIR}/05_freqk.sh" "${CONFIG_FILE}")
echo "Submitted 05_freqk.sh as job ${JOB5} (afterok:${JOB4})"

echo
echo "Pipeline submitted."
echo "  01_make_beds.sh : ${JOB1}"
echo "  02_run_hack.sh  : ${JOB2}"
echo "  03_run_shorts.sh: ${JOB3}"
echo "  04_make_vcf.sh  : ${JOB4}"
echo "  05_freqk.sh     : ${JOB5}"

