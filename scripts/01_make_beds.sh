#!/bin/bash
#SBATCH --job-name=make_beds
#SBATCH --output=logs/01_make_beds_%j.out
#SBATCH --error=logs/01_make_beds_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

# =============================================================================
# 01_make_beds.sh
# Purpose: Create VISOR HACk BED files for each SV size (DEL or INS)
#          Parameters come from a config file (default: config_sv_deletions.sh)
# BED format: chrom  start  end  svtype  info  stranddev
# =============================================================================

set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONFIG_FILE=${1:-"$(dirname "$0")/config_sv_deletions.sh"}
source "${CONFIG_FILE}"

mkdir -p "${BEDS}" logs

# ---------------------------------------------------------------------------
# Fast early-exit: if all BED files already exist AND encode the current
# SV_START_0, nothing to do.  A plain existence check is unsafe because
# re-using the same rep_id after a position change would silently propagate
# the old position into haplotypes, reads, VCFs and results.
# ---------------------------------------------------------------------------
_N_BED_DONE=0; _N_BED_TOTAL=0; _STALE=""
_check_bed_start() {
  local bed="$1" expected="$2"
  local got
  got=$(awk 'NR==1 {print $2; exit}' "${bed}" 2>/dev/null || true)
  [[ "${got}" == "${expected}" ]]
}
case "${SV_TYPE}" in
  "DEL") for _S in "${!DEL_SIZES[@]}"; do
           _N_BED_TOTAL=$((_N_BED_TOTAL+1))
           _BED="${BEDS}/hack_del_${_S}.bed"
           if [[ -s "${_BED}" ]]; then
             if _check_bed_start "${_BED}" "${SV_START_0}"; then
               _N_BED_DONE=$((_N_BED_DONE+1))
             else
               _STALE="${_STALE} ${_S}"
             fi
           fi
         done ;;
  "INS") for _S in "${!INS_SIZES[@]}"; do
           _N_BED_TOTAL=$((_N_BED_TOTAL+1))
           _BED="${BEDS}/hack_ins_${_S}.bed"
           if [[ -s "${_BED}" ]]; then
             if _check_bed_start "${_BED}" "${SV_START_0}"; then
               _N_BED_DONE=$((_N_BED_DONE+1))
             else
               _STALE="${_STALE} ${_S}"
             fi
           fi
         done ;;
esac
if [[ -n "${_STALE}" ]]; then
  echo "ERROR: existing BED file(s) for size(s)${_STALE} do not match current" >&2
  echo "       SV_START_0=${SV_START_0} in config.  Refusing to overwrite silently." >&2
  echo "       Delete the stale BEDs (under ${BEDS}) and all downstream artifacts" >&2
  echo "       (haplotypes_var, reads_var, vcf, results) for this rep, then re-run." >&2
  exit 1
fi
if [[ "${_N_BED_TOTAL}" -gt 0 && "${_N_BED_DONE}" -ge "${_N_BED_TOTAL}" ]]; then
  echo "[$(date)] Skipping 01_make_beds — all ${_N_BED_TOTAL} BED files already exist at POS=${SV_START_0} in ${BEDS}"
  exit 0
fi
echo "[$(date)] Found ${_N_BED_DONE}/${_N_BED_TOTAL} BED files at POS=${SV_START_0} — building missing ones."

source "$(mamba info --base)/etc/profile.d/conda.sh" && conda activate pang

case "${SV_TYPE}" in
  "DEL")
    echo "[$(date)] Creating deletion BED files in ${BEDS}"
    for SIZE in "${!DEL_SIZES[@]}"; do
      LEN=${DEL_SIZES[$SIZE]}
      END=$((SV_START_0 + LEN))
      BED="${BEDS}/hack_del_${SIZE}.bed"
      # Reuse only if the existing BED encodes the current SV_START_0
      if [[ -s "${BED}" ]] && _check_bed_start "${BED}" "${SV_START_0}"; then
        echo "[$(date)]   Reusing existing BED for ${SIZE}: ${BED}"
      else
        if [[ -s "${BED}" ]]; then
          echo "[$(date)]   Overwriting stale BED for ${SIZE}: ${BED}"
        else
          echo "[$(date)]   Writing BED for ${SIZE}: ${BED}"
        fi
        echo -e "${CHROM}\t${SV_START_0}\t${END}\tdeletion\tNone\t0" > "${BED}"
      fi
    done
    echo "[$(date)] Done. Deletion BED files created for active sizes:"
    for SIZE in "${!DEL_SIZES[@]}"; do
      ls -lh "${BEDS}/hack_del_${SIZE}.bed"
    done
    ;;
  "INS")
    echo "[$(date)] Creating insertion BED files in ${BEDS}"
    for SIZE in "${!INS_SIZES[@]}"; do
      LEN=${INS_SIZES[$SIZE]}
      BED="${BEDS}/hack_ins_${SIZE}.bed"
      if [[ -s "${BED}" ]]; then
        echo "[$(date)]   Reusing existing BED for ${SIZE}: ${BED}"
        continue
      fi
      # HACk expects column 5 to be the actual DNA sequence to insert (inf),
      # and column 6 to be the length of an additional random sequence at the breakpoint.
      # Use a reproducible, GC-aware random sequence (not poly-A) and set random seq length to 0.
      INS_SEQ=$(python3 - "$LEN" << 'PY'
import random, sys
bases = ['A','T','C','G']
weights = [0.32, 0.32, 0.18, 0.18]  # ~36% GC, reproducible
k = int(sys.argv[1])
print(''.join(random.choices(bases, weights=weights, k=k)))
PY
)
      # VISOR HACk uses s=e for insertions; start==end anchors a single base,
      # and the insertion is placed immediately after that anchor.
      START=${SV_START_0}
      END=${SV_START_0}
      echo -e "${CHROM}\t${START}\t${END}\tinsertion\t${INS_SEQ}\t0" > "${BED}"
    done
    echo "[$(date)] Done. Insertion BED files created for active sizes:"
    for SIZE in "${!INS_SIZES[@]}"; do
      ls -lh "${BEDS}/hack_ins_${SIZE}.bed"
    done
    ;;
  *)
    echo "01_make_beds.sh: unsupported SV_TYPE=${SV_TYPE} (expected DEL or INS)" >&2
    exit 1
    ;;
esac