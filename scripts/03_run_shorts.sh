#!/bin/bash -l
#SBATCH --job-name=visor_shorts_freq_all
#SBATCH --output=logs/03_run_shorts_%j.out
#SBATCH --error=logs/03_run_shorts_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
# =============================================================================
# 03_run_shorts.sh
# Purpose: Simulate pool-seq reads for freq benchmark using VISOR SHORtS
#          Parameters (freq, coverage, error, sizes, SV_TYPE) come from config
#          Can be run with or without sequencing errors via ERROR_RATE
# =============================================================================
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
source "$(mamba info --base)/etc/profile.d/conda.sh" && conda activate pang

CONFIG_FILE=${1:-"$(dirname "$0")/config_sv_deletions.sh"}
source "${CONFIG_FILE}"

# Derive labels for cov/freq/error so reads dirs are experiment-specific
COV_LABEL=$(awk -v c="${COVERAGE}" 'BEGIN{printf "cov%d", c}')
FREQ_LABEL=$(awk -v f="${FREQ}" 'BEGIN{printf "%.0f", f*100}')
# Encode error using digits after the decimal point in ERROR_RATE:
#   ERROR_RATE=0       -> err0
#   ERROR_RATE=0.001   -> err001
#   ERROR_RATE=0.01    -> err01
is_zero=$(awk -v e="${ERROR_RATE}" 'BEGIN{ if (e==0) print "yes"; else print "no"; }')
if [[ "${is_zero}" == "yes" ]]; then
  ERR_LABEL="0"
else
  if [[ "${ERROR_RATE}" == *.* ]]; then
    dec="${ERROR_RATE#*.}"
  else
    dec="${ERROR_RATE}"
  fi
  ERR_LABEL="${dec}"
fi
RUN_TAG="f${FREQ_LABEL}_err${ERR_LABEL}"

# WT clone (shared across all sizes; must exist — created by 02_run_hack)
CLONE_WT=${WT_CLONE}
if [[ ! -d "${CLONE_WT}" ]]; then
  echo "ERROR: WT clone directory missing: ${CLONE_WT}. Run 02_run_hack first (it creates _clone_WT)." >&2
  exit 1
fi
if ! compgen -G "${CLONE_WT}"/*.fa > /dev/null 2>&1; then
  echo "ERROR: No *.fa in WT clone: ${CLONE_WT}. Run 02_run_hack first." >&2
  exit 1
fi

mkdir -p logs "${BEDS}"

# Region BED (full Chr1, 5 columns)
CHR1_LEN=$(awk -v c="${CHROM}" '$1==c {print $2}' "${REF}.fai")
REGION_BED="${BEDS}/shorts_region_${CHROM}.bed"
echo -e "${CHROM}\t1\t${CHR1_LEN}\t100.0\t100.0" > "${REGION_BED}"
echo "[$(date)] Region BED: ${CHROM} 1-${CHR1_LEN}"

case "${SV_TYPE}" in
  "DEL")
    for SIZE in "${!DEL_SIZES[@]}"; do
      CLONE_DEL="${HAPS}/del_${SIZE}/HAP1"

      OUT="${READS}/${COV_LABEL}/freq_${SIZE}_${RUN_TAG}"
      rm -rf "${OUT}"
      mkdir -p "${OUT}"

      echo "[$(date)] Running VISOR SHORtS (DEL): 2 clones, freq=${FREQ}, size=${SIZE}, coverage=${COVERAGE}x, error=${ERROR_RATE}"

      VISOR SHORtS \
          -g "${REF}" \
          -s "${CLONE_DEL}" "${CLONE_WT}" \
          -b "${REGION_BED}" \
          -o "${OUT}" \
          --coverage "${COVERAGE}" \
          --clonefraction "$(awk -v f="${FREQ}" 'BEGIN{print f*100}')" "$(awk -v f="${FREQ}" 'BEGIN{print (1-f)*100}')" \
          --error "${ERROR_RATE}" \
          --fastq \
          --threads 8 || true

      echo "[$(date)] Done. Output in ${OUT}"
      ls -lh "${OUT}"
    done
    ;;
  "INS")
    for SIZE in "${!INS_SIZES[@]}"; do
      CLONE_INS="${HAPS}/ins_${SIZE}/HAP1"

      OUT="${READS}/${COV_LABEL}/freq_${SIZE}_${RUN_TAG}"
      rm -rf "${OUT}"
      mkdir -p "${OUT}"

      echo "[$(date)] Running VISOR SHORtS (INS): 2 clones, freq=${FREQ}, size=${SIZE}, coverage=${COVERAGE}x, error=${ERROR_RATE}"

      VISOR SHORtS \
          -g "${REF}" \
          -s "${CLONE_INS}" "${CLONE_WT}" \
          -b "${REGION_BED}" \
          -o "${OUT}" \
          --coverage "${COVERAGE}" \
          --clonefraction "$(awk -v f="${FREQ}" 'BEGIN{print f*100}')" "$(awk -v f="${FREQ}" 'BEGIN{print (1-f)*100}')" \
          --error "${ERROR_RATE}" \
          --fastq \
          --threads 8 || true

      echo "[$(date)] Done. Output in ${OUT}"
      ls -lh "${OUT}"
    done
    ;;
  *)
    echo "03_run_shorts.sh: unsupported SV_TYPE=${SV_TYPE} (expected DEL or INS)" >&2
    exit 1
    ;;
esac
