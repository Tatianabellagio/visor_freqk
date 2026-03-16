#!/bin/bash
# =============================================================================
# check_timing.sh
# Purpose: Report wall-clock time for each step of a submitted pipeline run.
#          Reads the .jobids file written by run_pipeline_var.sh (or any
#          run_pipeline*.sh that follows the same KEY=JOBID format).
#
# Usage:
#   bash scripts/check_timing.sh logs/pipeline_var_20260316_143200.jobids
#   bash scripts/check_timing.sh          # auto-selects latest .jobids file
# =============================================================================
set -euo pipefail

LOGS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../logs" && pwd)"

# ── Resolve jobids file ───────────────────────────────────────────────────────
if [[ $# -ge 1 ]]; then
  JOBIDS_FILE="$1"
else
  # Auto-select latest .jobids file in logs/
  JOBIDS_FILE=$(ls -t "${LOGS_DIR}"/pipeline_*.jobids 2>/dev/null | head -1 || true)
  if [[ -z "${JOBIDS_FILE}" ]]; then
    echo "No .jobids file found in ${LOGS_DIR}" >&2
    echo "Usage: bash $0 <path/to/pipeline_*.jobids>" >&2
    exit 1
  fi
  echo "Auto-selected: ${JOBIDS_FILE}"
fi

[[ -f "${JOBIDS_FILE}" ]] || { echo "File not found: ${JOBIDS_FILE}" >&2; exit 1; }

# ── Parse jobids file ─────────────────────────────────────────────────────────
declare -A STEP_TO_JOB
SUBMIT_TS=""
CONFIG=""

while IFS='=' read -r key val; do
  [[ -z "${key}" || "${key}" =~ ^# ]] && continue
  case "${key}" in
    SUBMIT_TS) SUBMIT_TS="${val}" ;;
    CONFIG)    CONFIG="${val}" ;;
    *)         STEP_TO_JOB["${key}"]="${val}" ;;
  esac
done < "${JOBIDS_FILE}"

# Ordered step list
STEPS=(
  00_prep_vcf
  00_apply_vcf
  01_make_beds
  02_run_hack_var
  03_run_shorts_var
  04_make_vcf_var
  05_freqk_var
)

# Collect all job IDs for a single sacct call
ALL_IDS=()
for step in "${STEPS[@]}"; do
  [[ -n "${STEP_TO_JOB[${step}]:-}" ]] && ALL_IDS+=("${STEP_TO_JOB[${step}]}")
done

if [[ ${#ALL_IDS[@]} -eq 0 ]]; then
  echo "No job IDs found in ${JOBIDS_FILE}" >&2
  exit 1
fi

IDS_CSV=$(IFS=,; echo "${ALL_IDS[*]}")

# ── Query sacct ───────────────────────────────────────────────────────────────
echo
echo "Pipeline run: ${SUBMIT_TS}"
echo "Config:       ${CONFIG}"
echo
printf "%-22s  %-10s  %-12s  %-10s  %-10s  %-8s\n" \
  "STEP" "JOB_ID" "STATE" "START" "END" "ELAPSED"
printf "%-22s  %-10s  %-12s  %-10s  %-10s  %-8s\n" \
  "$(printf '%.0s-' {1..22})" "----------" "------------" "----------" "----------" "--------"

# sacct returns one row per job step; use --parsable2 (| separator, no trailing |)
declare -A JOB_INFO
while IFS='|' read -r jobid jobname state elapsed start end _rest; do
  # skip sub-steps (jobid contains a dot, e.g. 7265.0)
  [[ "${jobid}" == *"."* ]] && continue
  JOB_INFO["${jobid}"]="${state}|${elapsed}|${start}|${end}"
done < <(sacct \
  --jobs="${IDS_CSV}" \
  --format="JobID,JobName,State,Elapsed,Start,End" \
  --noheader \
  --parsable2 2>/dev/null)

TOTAL_START=""
TOTAL_END=""

for step in "${STEPS[@]}"; do
  jid="${STEP_TO_JOB[${step}]:-}"
  if [[ -z "${jid}" ]]; then
    printf "%-22s  %-10s  %-12s  %-10s  %-10s  %-8s\n" \
      "${step}" "-" "-" "-" "-" "-"
    continue
  fi

  info="${JOB_INFO[${jid}]:-}"
  if [[ -z "${info}" ]]; then
    printf "%-22s  %-10s  %-12s  %-10s  %-10s  %-8s\n" \
      "${step}" "${jid}" "PENDING/UNK" "-" "-" "-"
    continue
  fi

  IFS='|' read -r state elapsed start end <<< "${info}"

  # Track pipeline-level start/end
  if [[ "${start}" != "None" && "${start}" != "Unknown" ]]; then
    [[ -z "${TOTAL_START}" || "${start}" < "${TOTAL_START}" ]] && TOTAL_START="${start}"
  fi
  if [[ "${end}" != "None" && "${end}" != "Unknown" ]]; then
    [[ -z "${TOTAL_END}"   || "${end}"   > "${TOTAL_END}"   ]] && TOTAL_END="${end}"
  fi

  # Shorten timestamps to HH:MM:SS for readability
  start_short=$(echo "${start}" | grep -oE '[0-9]{2}:[0-9]{2}:[0-9]{2}' || echo "${start}")
  end_short=$(  echo "${end}"   | grep -oE '[0-9]{2}:[0-9]{2}:[0-9]{2}' || echo "${end}")

  printf "%-22s  %-10s  %-12s  %-10s  %-10s  %-8s\n" \
    "${step}" "${jid}" "${state}" "${start_short}" "${end_short}" "${elapsed}"
done

echo
# Total wall time (start of first job → end of last job)
if [[ -n "${TOTAL_START}" && -n "${TOTAL_END}" ]]; then
  WALL=$(python3 -c "
from datetime import datetime
fmt = '%Y-%m-%dT%H:%M:%S'
s = datetime.strptime('${TOTAL_START}', fmt)
e = datetime.strptime('${TOTAL_END}',   fmt)
d = e - s
h, rem = divmod(int(d.total_seconds()), 3600)
m, sec = divmod(rem, 60)
print(f'{h:02d}:{m:02d}:{sec:02d}')
" 2>/dev/null || echo "?")
  echo "Total wall time (first start → last end): ${WALL}"
fi
echo "Raw sacct for full detail:"
echo "  sacct -j ${IDS_CSV} --format=JobID,JobName,State,Elapsed,Start,End,MaxRSS"
