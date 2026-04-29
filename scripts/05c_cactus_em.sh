#!/bin/bash
#SBATCH --job-name=cactus_em
#SBATCH --output=logs/05c_cactus_em_%j.out
#SBATCH --error=logs/05c_cactus_em_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

# =============================================================================
# 05c_cactus_em.sh
# Replaces 05_freqk_var.sh for the cactus_em pipeline. Runs our cactus
# pangenome k-mer EM on the SHORtS-pooled reads, in both global and
# per-window modes for comparison.
#
# Output: ${RESULTS}/cov${COV}/var_del_${SIZE}_<tag>/<mode>.tsv
#         (per-VCF-record alt_freq table, schema matching freqk's output)
# =============================================================================

set -euo pipefail
CONFIG_FILE=${1:-"$(dirname "$0")/config_cactus_em.sh"}
source "${CONFIG_FILE}"

mkdir -p logs ${RESULTS}

FREQ_LABEL=$(python3 -c "print(round(${SV_FREQ}*100))")
ERR_LABEL=$(python3 -c "
e='${ERROR_RATE}'
print('0' if float(e)==0 else (e.split('.')[-1] if '.' in e else e))")
RUN_TAG="n${N_SAMPLES}_f${FREQ_LABEL}_err${ERR_LABEL}"

for SIZE in "${!DEL_SIZES[@]}"; do
    READS_DIR=${READS_VAR}/cov${COVERAGE}/var_del_${SIZE}_${RUN_TAG}
    OUT_DIR=${RESULTS}/cov${COVERAGE}/var_del_${SIZE}_${RUN_TAG}
    mkdir -p ${OUT_DIR}

    if [ ! -s ${READS_DIR}/r1.fq ] || [ ! -s ${READS_DIR}/r2.fq ]; then
        echo "ERROR: reads missing in ${READS_DIR}" >&2; exit 1
    fi

    for MODE in global window; do
        OUT=${OUT_DIR}/cactus_em_${MODE}.tsv
        if [ -s ${OUT} ]; then echo "[$(date)] $MODE exists: ${OUT}"; continue; fi
        echo ""
        echo "============================================================"
        echo "[$(date)] cactus_em ${SIZE} ${MODE}"
        echo "============================================================"
        ${PYTHON_HAPFM} -u ${POOLFREQ}/src/per_sample_driver.py \
            --cn-kmer-prefix ${CN_KMER_PREFIX} \
            --cn-var ${CN_VAR} \
            --cn-var-meta ${CN_VAR_META} \
            --reads ${READS_DIR}/r1.fq ${READS_DIR}/r2.fq \
            --sample ${POS_LABEL}_${SIZE}_${MODE} \
            --out ${OUT} \
            --threads 8 \
            --block-mode ${MODE} \
            --window-bp 200000
    done
done

echo "[$(date)] DONE. Results:"
find ${RESULTS} -name "cactus_em_*.tsv" -ls
