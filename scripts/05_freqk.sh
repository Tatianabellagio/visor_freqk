#!/bin/bash -l
#SBATCH --job-name=freqk_all
#SBATCH --output=logs/05_freqk_%j.out
#SBATCH --error=logs/05_freqk_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
# =============================================================================
# 05_freqk.sh
# Purpose: Full freqk pipeline for each SV-size VCF:
#          index + var-dedup + ref-dedup + count + call
#          Parameters (paths, k, freq/error tags, SV_TYPE) come from config file
# =============================================================================
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
eval "$(conda shell.bash hook)"
conda activate freqk_build

CONFIG_FILE=${1:-"$(dirname "$0")/config_sv_deletions.sh"}
source "${CONFIG_FILE}"

mkdir -p "${RESULTS}"

# Derive the same RUN_TAG as in 03_run_shorts.sh (cov/freq/error)
FREQ_LABEL=$(awk -v f="${FREQ}" 'BEGIN{printf "%.0f", f*100}')
is_zero_rt=$(awk -v e="${ERROR_RATE}" 'BEGIN{ if (e==0) print "yes"; else print "no"; }')
if [[ "${is_zero_rt}" == "yes" ]]; then
  ERR_LABEL="0"
else
  if [[ "${ERROR_RATE}" == *.* ]]; then
    dec_rt="${ERROR_RATE#*.}"
  else
    dec_rt="${ERROR_RATE}"
  fi
  ERR_LABEL="${dec_rt}"
fi
RUN_TAG="f${FREQ_LABEL}_err${ERR_LABEL}"

# Labels for requested results folder structure:
# results/coverage_sequencingerror/sv_type/sv_length/sv_freq/freqk_k/
COV_LABEL=$(awk -v c="${COVERAGE}" 'BEGIN{printf "cov%d", c}')
# Encode full error rate using the digits after the decimal point in ERROR_RATE:
#   ERROR_RATE=0       -> err0
#   ERROR_RATE=0.001   -> err001
#   ERROR_RATE=0.01    -> err01
#   ERROR_RATE=0.123   -> err123
is_zero=$(awk -v e="${ERROR_RATE}" 'BEGIN{ if (e==0) print "yes"; else print "no"; }')
if [[ "${is_zero}" == "yes" ]]; then
  ERR_LABEL_DIR="err0"
else
  if [[ "${ERROR_RATE}" == *.* ]]; then
    dec="${ERROR_RATE#*.}"
  else
    dec="${ERROR_RATE}"
  fi
  ERR_LABEL_DIR="err${dec}"
fi
FREQ_LABEL_DIR=$(awk -v f="${FREQ}" 'BEGIN{printf "f%.0f", f*100}')
K_LABEL="k${K}"

case "${SV_TYPE}" in
  "DEL")
    for SIZE in "${!DEL_SIZES[@]}"; do
      VCF=${VCF_DIR}/del_${SIZE}.vcf.gz

      # results/coverage_sequencingerror/sv_type/sv_length/sv_freq/freqk_k/
      BASE_DIR=${RESULTS}/${COV_LABEL}_${ERR_LABEL_DIR}/${SV_TYPE}/${SIZE}/${FREQ_LABEL_DIR}/${K_LABEL}

      INDEX=${BASE_DIR}/del_${SIZE}.k${K}.freqk.index
      VAR_INDEX=${BASE_DIR}/del_${SIZE}.k${K}.freqk.var_index
      REF_INDEX=${BASE_DIR}/del_${SIZE}.k${K}.freqk.ref_index

      READS_DIR=${READS}/${COV_LABEL}/freq_${SIZE}_${RUN_TAG}
      RESULTS_DIR=${BASE_DIR}
      mkdir -p "${RESULTS_DIR}"

      READS_COMBINED=${READS_DIR}/all.fq
      COUNTS_BY_ALLELE=${RESULTS_DIR}/freq_${SIZE}_${RUN_TAG}.counts_by_allele.k${K}.tsv
      RAW_COUNTS=${RESULTS_DIR}/freq_${SIZE}_${RUN_TAG}.raw_kmer_counts.k${K}.tsv
      AF_OUT=${RESULTS_DIR}/freq_${SIZE}_${RUN_TAG}.allele_frequencies.k${K}.tsv

      echo
      echo "==== Size ${SIZE} (RUN_TAG=${RUN_TAG}, SV_TYPE=DEL) ===="

      [[ -s "$REF" && -s "${REF}.fai" ]] || { echo "Missing FASTA or fai" >&2; exit 1; }
      [[ -s "$VCF"   && -s "${VCF}.tbi"   ]] || { echo "Missing VCF or tbi for ${SIZE}"  >&2; continue; }
      [[ -s "$READS_DIR/r1.fq" && -s "$READS_DIR/r2.fq" ]] || { echo "Missing reads for ${SIZE} at ${READS_DIR}" >&2; continue; }

      echo "[$(date)] Building freqk index (k=${K}) for ${SIZE}"
      "$FREQK" index --fasta "$REF" --vcf "$VCF" --output "$INDEX" --kmer "$K"
      ls -lh "$INDEX"

      echo "[$(date)] Combining r1.fq + r2.fq -> all.fq (${SIZE})"
      cat "${READS_DIR}/r1.fq" "${READS_DIR}/r2.fq" > "$READS_COMBINED"

      step() {
        local label="$1"; shift
        echo
        echo "[$(date)] === $label (${SIZE}) ==="
        echo "+ $*"
        "$@"
      }

      step "var-dedup" \
        "$FREQK" var-dedup --index "$INDEX" --output "$VAR_INDEX"
      ls -lh "$VAR_INDEX"

      step "ref-dedup" \
        "$FREQK" ref-dedup -i "$VAR_INDEX" -o "$REF_INDEX" -f "$REF" --vcf "$VCF"
      ls -lh "$REF_INDEX"

      step "count" \
        "$FREQK" count \
          --index "$REF_INDEX" \
          --reads "$READS_COMBINED" \
          --nthreads "${SLURM_CPUS_PER_TASK}" \
          --freq-output "$COUNTS_BY_ALLELE" \
          --count-output "$RAW_COUNTS"

      step "call" \
        "$FREQK" call \
          --index "$REF_INDEX" \
          --counts "$COUNTS_BY_ALLELE" \
          --output "$AF_OUT"

      echo "[$(date)] Done for ${SIZE}. AF:"
      cat "$AF_OUT"
    done
    ;;
  "INS")
    for SIZE in "${!INS_SIZES[@]}"; do
      VCF=${VCF_DIR}/ins_${SIZE}.vcf.gz

      BASE_DIR=${RESULTS}/${COV_LABEL}_${ERR_LABEL_DIR}/${SV_TYPE}/${SIZE}/${FREQ_LABEL_DIR}/${K_LABEL}

      INDEX=${BASE_DIR}/ins_${SIZE}.k${K}.freqk.index
      VAR_INDEX=${BASE_DIR}/ins_${SIZE}.k${K}.freqk.var_index
      REF_INDEX=${BASE_DIR}/ins_${SIZE}.k${K}.freqk.ref_index

      READS_DIR=${READS}/${COV_LABEL}/freq_${SIZE}_${RUN_TAG}
      RESULTS_DIR=${BASE_DIR}
      mkdir -p "${RESULTS_DIR}"

      READS_COMBINED=${READS_DIR}/all.fq
      COUNTS_BY_ALLELE=${RESULTS_DIR}/freq_${SIZE}_${RUN_TAG}.counts_by_allele.k${K}.tsv
      RAW_COUNTS=${RESULTS_DIR}/freq_${SIZE}_${RUN_TAG}.raw_kmer_counts.k${K}.tsv
      AF_OUT=${RESULTS_DIR}/freq_${SIZE}_${RUN_TAG}.allele_frequencies.k${K}.tsv

      echo
      echo "==== Size ${SIZE} (RUN_TAG=${RUN_TAG}, SV_TYPE=INS) ===="

      [[ -s "$REF" && -s "${REF}.fai" ]] || { echo "Missing FASTA or fai" >&2; exit 1; }
      [[ -s "$VCF"   && -s "${VCF}.tbi"   ]] || { echo "Missing VCF or tbi for ${SIZE}"  >&2; continue; }
      [[ -s "$READS_DIR/r1.fq" && -s "$READS_DIR/r2.fq" ]] || { echo "Missing reads for ${SIZE} at ${READS_DIR}" >&2; continue; }

      echo "[$(date)] Building freqk index (k=${K}) for ${SIZE}"
      "$FREQK" index --fasta "$REF" --vcf "$VCF" --output "$INDEX" --kmer "$K"
      ls -lh "$INDEX"

      echo "[$(date)] Combining r1.fq + r2.fq -> all.fq (${SIZE})"
      cat "${READS_DIR}/r1.fq" "${READS_DIR}/r2.fq" > "$READS_COMBINED"

      step() {
        local label="$1"; shift
        echo
        echo "[$(date)] === $label (${SIZE}) ==="
        echo "+ $*"
        "$@"
      }

      step "var-dedup" \
        "$FREQK" var-dedup --index "$INDEX" --output "$VAR_INDEX"
      ls -lh "$VAR_INDEX"

      step "ref-dedup" \
        "$FREQK" ref-dedup -i "$VAR_INDEX" -o "$REF_INDEX" -f "$REF" --vcf "$VCF"
      ls -lh "$REF_INDEX"

      step "count" \
        "$FREQK" count \
          --index "$REF_INDEX" \
          --reads "$READS_COMBINED" \
          --nthreads "${SLURM_CPUS_PER_TASK}" \
          --freq-output "$COUNTS_BY_ALLELE" \
          --count-output "$RAW_COUNTS"

      step "call" \
        "$FREQK" call \
          --index "$REF_INDEX" \
          --counts "$COUNTS_BY_ALLELE" \
          --output "$AF_OUT"

      echo "[$(date)] Done for ${SIZE}. AF:"
      cat "$AF_OUT"
    done
    ;;
  *)
    echo "05_freqk.sh: unsupported SV_TYPE=${SV_TYPE} (expected DEL or INS)" >&2
    exit 1
    ;;
esac