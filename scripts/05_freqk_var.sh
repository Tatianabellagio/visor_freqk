#!/bin/bash -l
#SBATCH --job-name=freqk_var
#SBATCH --output=logs/05_freqk_var_%j.out
#SBATCH --error=logs/05_freqk_var_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
# =============================================================================
# 05_freqk_var.sh
# Purpose: Full freqk pipeline on reads produced by 03_run_shorts_var.sh.
#          Matches the var-pipeline directory naming:
#            reads_var/del/<pos>/cov<X>/var_del_<SIZE>_n<N>_f<F>_err<E>/
#          rather than the original pipeline's freq_<SIZE>_f<F>_err<E>/ layout.
#
# Usage:
#   sbatch scripts/05_freqk_var.sh [config_file]
#   (default config: config_sv_var_deletions.sh)
# =============================================================================
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONFIG_FILE=${1:-"$(dirname "$0")/config_sv_var_deletions.sh"}
source "${CONFIG_FILE}"

mkdir -p "${RESULTS}"

# ---------------------------------------------------------------------------
# Fast early-exit: if all allele-frequency files already exist, nothing to do.
# Runs before conda activation to avoid unnecessary overhead.
# ---------------------------------------------------------------------------
_FL=$(python3 -c "print(round(${SV_FREQ}*100))")
_EL=$(python3 -c "e='${ERROR_RATE}'; print('0' if float(e)==0 else (e.split('.')[-1] if '.' in e else e))")
_RTAG="n${N_SAMPLES}_f${_FL}_err${_EL}"
_N_AF_DONE=0; _N_AF_TOTAL=0
for _SZ in "${!DEL_SIZES[@]}"; do
  _N_AF_TOTAL=$((_N_AF_TOTAL + 1))
  _BDIR="${RESULTS}/cov${COVERAGE}_err${_EL}/${_SZ}/n${N_SAMPLES}/f${_FL}/k${K}"
  _AF="${_BDIR}/var_del_${_SZ}_${_RTAG}.allele_frequencies.k${K}.tsv"
  [[ -s "${_AF}" ]] && _N_AF_DONE=$((_N_AF_DONE + 1))
done
if [[ "${_N_AF_TOTAL}" -gt 0 && "${_N_AF_DONE}" -ge "${_N_AF_TOTAL}" ]]; then
  echo "[$(date)] Skipping 05_freqk_var — all ${_N_AF_TOTAL} allele-frequency files already exist"
  exit 0
fi
echo "[$(date)] Found ${_N_AF_DONE}/${_N_AF_TOTAL} allele-frequency files — running missing ones."

eval "$(conda shell.bash hook)"
conda activate freqk_build

# ---------------------------------------------------------------------------
# Build label components  (must mirror 03_run_shorts_var.sh exactly)
# ---------------------------------------------------------------------------
FREQ_LABEL=$(python3 -c "print(round(${SV_FREQ}*100))")
ERR_LABEL=$(python3 -c "
e='${ERROR_RATE}'
if float(e)==0: print('0')
else: print(e.split('.')[-1] if '.' in e else e)
")
COV_LABEL="cov${COVERAGE}"

# Run-tag used by 03_run_shorts_var.sh: n<N>_f<F>_err<E>
RUN_TAG="n${N_SAMPLES}_f${FREQ_LABEL}_err${ERR_LABEL}"

# Directory-level labels for results folder structure
ERR_LABEL_DIR="err${ERR_LABEL}"
FREQ_LABEL_DIR="f${FREQ_LABEL}"
K_LABEL="k${K}"

echo "[$(date)] Config:   ${CONFIG_FILE}"
echo "[$(date)] RUN_TAG:  ${RUN_TAG}"
echo "[$(date)] RESULTS:  ${RESULTS}"

# ---------------------------------------------------------------------------
case "${SV_TYPE}" in
  "DEL")
    for SIZE in "${!DEL_SIZES[@]}"; do

      VCF=${VCF_DIR}/del_${SIZE}.vcf.gz

      # results/cov<X>_err<E>/<SIZE>/n<N>/f<F>/k<K>/
      N_LABEL="n${N_SAMPLES}"
      BASE_DIR=${RESULTS}/${COV_LABEL}_${ERR_LABEL_DIR}/${SIZE}/${N_LABEL}/${FREQ_LABEL_DIR}/${K_LABEL}

      INDEX=${BASE_DIR}/del_${SIZE}.k${K}.freqk.index
      VAR_INDEX=${BASE_DIR}/del_${SIZE}.k${K}.freqk.var_index
      REF_INDEX=${BASE_DIR}/del_${SIZE}.k${K}.freqk.ref_index

      # Path written by 03_run_shorts_var.sh
      READS_DIR=${READS}/${COV_LABEL}/var_del_${SIZE}_${RUN_TAG}
      mkdir -p "${BASE_DIR}"

      READS_COMBINED=${READS_DIR}/all.fq
      COUNTS_BY_ALLELE=${BASE_DIR}/var_del_${SIZE}_${RUN_TAG}.counts_by_allele.k${K}.tsv
      RAW_COUNTS=${BASE_DIR}/var_del_${SIZE}_${RUN_TAG}.raw_kmer_counts.k${K}.tsv
      AF_OUT=${BASE_DIR}/var_del_${SIZE}_${RUN_TAG}.allele_frequencies.k${K}.tsv
      RUNTIMES=${BASE_DIR}/var_del_${SIZE}_${RUN_TAG}.runtimes.k${K}.tsv

      echo
      echo "==== DEL ${SIZE}  (RUN_TAG=${RUN_TAG}) ===="
      echo "     VCF:      ${VCF}"
      echo "     READS:    ${READS_DIR}"
      echo "     RESULTS:  ${BASE_DIR}"

      if [[ -s "${AF_OUT}" ]]; then
        echo "[$(date)] Skipping DEL ${SIZE} — allele-frequency file already exists: ${AF_OUT}"
        continue
      fi

      # Preflight checks
      [[ -s "${REF}" && -s "${REF}.fai" ]] \
        || { echo "ERROR: Missing reference FASTA or .fai: ${REF}" >&2; exit 1; }
      [[ -s "${VCF}" && -s "${VCF}.tbi" ]] \
        || { echo "SKIP: Missing VCF or .tbi for size ${SIZE}: ${VCF}" >&2; continue; }
      [[ -s "${READS_DIR}/r1.fq" && -s "${READS_DIR}/r2.fq" ]] \
        || { echo "ERROR: Missing reads at ${READS_DIR}" >&2
             echo "       Expected r1.fq and r2.fq — did 03_run_shorts_var.sh succeed?" >&2
             exit 1; }

      # Per-step runtime TSV (sibling of AF output): method size step elapsed_s cached
      printf "method\tsize\tstep\telapsed_s\tcached\n" > "${RUNTIMES}.tmp"
      mv -f "${RUNTIMES}.tmp" "${RUNTIMES}"

      step() {
        local label="$1"; shift
        local _t0=${SECONDS}
        echo
        echo "[$(date)] === ${label} (DEL ${SIZE}) ==="
        echo "+ $*"
        "$@"
        local _el=$(( SECONDS - _t0 ))
        printf "freqk\t%s\t%s\t%d\t0\n" "${SIZE}" "${label}" "${_el}" >> "${RUNTIMES}"
      }

      # --- index ---
      step "index" \
        "${FREQK}" index --fasta "${REF}" --vcf "${VCF}" --output "${INDEX}" --kmer "${K}"
      ls -lh "${INDEX}"

      # --- combine reads ---
      # Write to a temp file then atomic mv so concurrent k-jobs sharing the
      # same reads dir don't truncate the file while another job is reading it.
      if [[ ! -s "${READS_COMBINED}" ]]; then
        echo "[$(date)] Combining r1.fq + r2.fq → all.fq (${SIZE})"
        _t0=${SECONDS}
        cat "${READS_DIR}/r1.fq" "${READS_DIR}/r2.fq" > "${READS_COMBINED}.tmp"
        mv -f "${READS_COMBINED}.tmp" "${READS_COMBINED}"
        printf "freqk\t%s\tcombine_reads\t%d\t0\n" "${SIZE}" $(( SECONDS - _t0 )) >> "${RUNTIMES}"
      else
        echo "[$(date)] Reusing existing all.fq (${SIZE})"
        printf "freqk\t%s\tcombine_reads\t0\t1\n" "${SIZE}" >> "${RUNTIMES}"
      fi

      step "var-dedup" \
        "${FREQK}" var-dedup --index "${INDEX}" --output "${VAR_INDEX}"
      ls -lh "${VAR_INDEX}"

      step "ref-dedup" \
        "${FREQK}" ref-dedup -i "${VAR_INDEX}" -o "${REF_INDEX}" -f "${REF}" --vcf "${VCF}"
      ls -lh "${REF_INDEX}"

      step "count" \
        "${FREQK}" count \
          --index "${REF_INDEX}" \
          --reads "${READS_COMBINED}" \
          --nthreads "${SLURM_CPUS_PER_TASK}" \
          --freq-output "${COUNTS_BY_ALLELE}" \
          --count-output "${RAW_COUNTS}"

      step "call" \
        "${FREQK}" call \
          --index "${REF_INDEX}" \
          --counts "${COUNTS_BY_ALLELE}" \
          --output "${AF_OUT}"

      echo
      echo "[$(date)] Done for DEL ${SIZE}. Allele frequencies:"
      cat "${AF_OUT}"
    done
    ;;
  *)
    echo "05_freqk_var.sh: unsupported SV_TYPE=${SV_TYPE} (expected DEL)" >&2
    exit 1
    ;;
esac
