#!/bin/bash
#SBATCH --job-name=cactus_shorts
#SBATCH --output=logs/03_run_shorts_cactus_%j.out
#SBATCH --error=logs/03_run_shorts_cactus_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# =============================================================================
# 03_run_shorts_cactus.sh
# Run VISOR SHORtS to pool reads from N cactus founders, where the first
# N_SV carry the deletion. Sample list comes from samples.txt.
# =============================================================================

set -euo pipefail
CONFIG_FILE=${1:-"$(dirname "$0")/config_cactus_em.sh"}
source "${CONFIG_FILE}"

mkdir -p logs ${BEDS}
source "$(mamba info --base)/etc/profile.d/conda.sh" && conda activate pang

# Sample list
SAMPLES_FILE=$(dirname "${HAPS_WT}")/samples.txt
mapfile -t SAMPLES < ${SAMPLES_FILE}
N_SV=$(python3 -c "print(round(${N_SAMPLES} * ${SV_FREQ}))")
N_WT=$(( N_SAMPLES - N_SV ))

# Helper: equal fractions summing to exactly 100
make_fracs() {
    python3 - "$1" << 'PY'
import sys
n = int(sys.argv[1])
base = round(100.0 / n, 4)
fracs = [base] * n
fracs[-1] = round(100.0 - sum(fracs[:-1]), 4)
print('\n'.join(str(x) for x in fracs))
PY
}

for SIZE in "${!DEL_SIZES[@]}"; do
    # Build clone dir list: SV first, then WT (matches existing visor_freqk order)
    CLONE_DIRS=()
    for SAMPLE in "${SAMPLES[@]:0:${N_SV}}"; do
        CLONE_DIRS+=("${HAPS_SV}/s_${SAMPLE}_sv_del_${SIZE}")
    done
    for SAMPLE in "${SAMPLES[@]:${N_SV}}"; do
        CLONE_DIRS+=("${HAPS_WT}/s_${SAMPLE}")
    done
    N_CLONES=${#CLONE_DIRS[@]}
    mapfile -t FRACS < <(make_fracs ${N_CLONES})
    echo "[$(date)] N=${N_CLONES} clones, fractions: ${FRACS[*]}"

    FREQ_LABEL=$(python3 -c "print(round(${SV_FREQ}*100))")
    ERR_LABEL=$(python3 -c "
e='${ERROR_RATE}'
print('0' if float(e)==0 else (e.split('.')[-1] if '.' in e else e))")
    RUN_TAG="n${N_SAMPLES}_f${FREQ_LABEL}_err${ERR_LABEL}"
    OUT="${READS_VAR}/cov${COVERAGE}/var_del_${SIZE}_${RUN_TAG}"

    if [ -s ${OUT}/r1.fq ] && [ -s ${OUT}/r2.fq ]; then
        echo "[$(date)] Reusing existing reads in ${OUT}"
        continue
    fi
    rm -rf ${OUT}; mkdir -p ${OUT}

    echo "[$(date)] VISOR SHORtS: DEL ${SIZE}, cov=${COVERAGE}× on ${CHROM}"
    VISOR SHORtS \
        -g ${REF} \
        -s ${CLONE_DIRS[@]} \
        -b ${REGION_BED_PATH} \
        -o ${OUT} \
        --coverage ${COVERAGE} \
        --clonefraction ${FRACS[@]} \
        --error ${ERROR_RATE} \
        --length ${READ_LEN} \
        --fastq \
        --threads 8 || true

    if [ ! -s ${OUT}/r1.fq ] || [ ! -s ${OUT}/r2.fq ]; then
        echo "ERROR: SHORtS missing r1.fq/r2.fq in ${OUT}" >&2; exit 1
    fi
    echo "[$(date)] Done: ${OUT}"
    ls -lh ${OUT}
done
