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
source "$(mamba info --base)/etc/profile.d/conda.sh" && conda activate pang

CONFIG_FILE=${1:-"$(dirname "$0")/config_sv_deletions.sh"}
source "${CONFIG_FILE}"

mkdir -p "${BEDS}" logs

case "${SV_TYPE}" in
  "DEL")
    echo "[$(date)] Creating deletion BED files in ${BEDS}"
    for SIZE in "${!DEL_SIZES[@]}"; do
      LEN=${DEL_SIZES[$SIZE]}
      END=$((SV_START_0 + LEN))
      BED="${BEDS}/hack_del_${SIZE}.bed"
      # If a BED with these parameters already exists, reuse it
      if [[ -s "${BED}" ]]; then
        echo "[$(date)]   Reusing existing BED for ${SIZE}: ${BED}"
      else
        echo "[$(date)]   Writing BED for ${SIZE}: ${BED}"
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