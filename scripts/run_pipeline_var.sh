#!/bin/bash
# =============================================================================
# run_pipeline_var.sh
# Purpose: Submit the variation-aware pipeline as a chained series of Slurm jobs.
#
# Job chain:
#   00_prep_vcf ──┬── 00_apply_vcf ──┐
#                 └── 01_make_beds ──┴── 02_run_hack_var ── 03_run_shorts_var ──┐
#   04_make_vcf_var (no deps, runs immediately) ──────────────────────────────────┴── 05_freqk_var
#
# Usage:
#   cd /home/tbellagio/scratch/visor_freqk
#   bash scripts/run_pipeline_var.sh [config_file]
#   (default config: scripts/config_sv_var_deletions.sh)
# =============================================================================
set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE=${1:-"${SCRIPTS_DIR}/config_sv_var_deletions.sh"}

if [[ ! -f "${CONFIG_FILE}" ]]; then
  echo "Config file not found: ${CONFIG_FILE}" >&2
  exit 1
fi

SUBMIT_TS=$(date +%Y%m%d_%H%M%S)
JOBIDS_FILE="${SCRIPTS_DIR}/../logs/pipeline_var_${SUBMIT_TS}.jobids"
mkdir -p "$(dirname "${JOBIDS_FILE}")"

echo "Using config: ${CONFIG_FILE}"
echo "Submit time:  ${SUBMIT_TS}"

# Step 00_prep: rename VCF chromosomes once ("1" → "Chr1")
JOB_PREP=$(sbatch --parsable "${SCRIPTS_DIR}/00_prep_vcf.sh" "${CONFIG_FILE}")
echo "Submitted 00_prep_vcf.sh      as job ${JOB_PREP}"

# Step 00_apply + 01 both depend on prep (run in parallel with each other)
JOB0=$(sbatch --parsable --dependency=afterok:${JOB_PREP} "${SCRIPTS_DIR}/00_apply_vcf.sh" "${CONFIG_FILE}")
echo "Submitted 00_apply_vcf.sh     as job ${JOB0} (afterok:${JOB_PREP})"

JOB1=$(sbatch --parsable --dependency=afterok:${JOB_PREP} "${SCRIPTS_DIR}/01_make_beds.sh" "${CONFIG_FILE}")
echo "Submitted 01_make_beds.sh     as job ${JOB1} (afterok:${JOB_PREP})"

# Step 02var: VISOR HACk on SV-bearing samples (needs 00_apply + 01)
JOB2=$(sbatch --parsable --dependency=afterok:${JOB0}:${JOB1} "${SCRIPTS_DIR}/02_run_hack_var.sh" "${CONFIG_FILE}")
echo "Submitted 02_run_hack_var.sh  as job ${JOB2} (afterok:${JOB0}:${JOB1})"

# Step 03var: SHORtS pool all N clones (needs 02var)
JOB3=$(sbatch --parsable --dependency=afterok:${JOB2} "${SCRIPTS_DIR}/03_run_shorts_var.sh" "${CONFIG_FILE}")
echo "Submitted 03_run_shorts_var.sh as job ${JOB3} (afterok:${JOB2})"

# Step 04var: build SV VCFs — no dependency on haplotypes or reads, runs immediately
JOB4=$(sbatch --parsable "${SCRIPTS_DIR}/04_make_vcf_var.sh" "${CONFIG_FILE}")
echo "Submitted 04_make_vcf_var.sh   as job ${JOB4} (no dependency)"

# Step 05var: freqk — needs both reads (JOB3) and VCF (JOB4)
JOB5=$(sbatch --parsable --dependency=afterok:${JOB3}:${JOB4} "${SCRIPTS_DIR}/05_freqk_var.sh" "${CONFIG_FILE}")
echo "Submitted 05_freqk_var.sh      as job ${JOB5} (afterok:${JOB3}:${JOB4})"

# Save job IDs for timing queries
cat > "${JOBIDS_FILE}" << EOF
SUBMIT_TS=${SUBMIT_TS}
CONFIG=${CONFIG_FILE}
00_prep_vcf=${JOB_PREP}
00_apply_vcf=${JOB0}
01_make_beds=${JOB1}
02_run_hack_var=${JOB2}
03_run_shorts_var=${JOB3}
04_make_vcf_var=${JOB4}
05_freqk_var=${JOB5}
EOF

echo
echo "Pipeline submitted."
echo "  00_prep_vcf.sh        : ${JOB_PREP}"
echo "  00_apply_vcf.sh       : ${JOB0}"
echo "  01_make_beds.sh       : ${JOB1}"
echo "  02_run_hack_var.sh    : ${JOB2}"
echo "  03_run_shorts_var.sh  : ${JOB3}"
echo "  04_make_vcf_var.sh    : ${JOB4}"
echo "  05_freqk_var.sh       : ${JOB5}"
echo
echo "Job IDs saved to: ${JOBIDS_FILE}"
echo "Check timing once done:"
echo "  bash ${SCRIPTS_DIR}/check_timing.sh ${JOBIDS_FILE}"
