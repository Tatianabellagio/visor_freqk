#!/bin/bash
#SBATCH --job-name=cactus_hack
#SBATCH --output=logs/02_run_hack_cactus_%j.out
#SBATCH --error=logs/02_run_hack_cactus_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# =============================================================================
# 02_run_hack_cactus.sh
# Same logic as 02_run_hack_var.sh but reads SAMPLES from samples.txt
# (written by 00_setup_cactus_clones.sh) instead of bcftools query of a VCF.
# =============================================================================

set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
CONFIG_FILE=${1:-"$(dirname "$0")/config_cactus_em.sh"}
source "${CONFIG_FILE}"

mkdir -p logs "${HAPS_SV}"
source "$(mamba info --base)/etc/profile.d/conda.sh" && conda activate pang

# Read sample list (from 00_setup_cactus_clones.sh)
SAMPLES_FILE=$(dirname "${HAPS_WT}")/samples.txt
mapfile -t SAMPLES < ${SAMPLES_FILE}
N_SV=$(python3 -c "print(round(${N_SAMPLES} * ${SV_FREQ}))")
N_WT=$(( N_SAMPLES - N_SV ))
echo "[$(date)] Samples: ${SAMPLES[*]}"
echo "[$(date)] N_SAMPLES=${N_SAMPLES}  N_SV=${N_SV}  N_WT=${N_WT}"

# SV-bearing samples are the FIRST N_SV from the sorted list
SV_SAMPLES=("${SAMPLES[@]:0:${N_SV}}")
echo "[$(date)] SV carriers: ${SV_SAMPLES[*]}"

# Loop over deletion sizes
for SIZE in "${!DEL_SIZES[@]}"; do
    LEN=${DEL_SIZES[$SIZE]}
    BED=${BEDS}/del_${SIZE}.bed
    mkdir -p ${BEDS}
    SV_END_0=$((SV_START_0 + LEN))
    echo -e "${CHROM}\t${SV_START_0}\t${SV_END_0}\tdeletion\tNone\t0" > ${BED}
    echo "[$(date)] HACk BED for ${SIZE}: ${BED}"

    for SAMPLE in "${SV_SAMPLES[@]}"; do
        SAMPLE_FA=${HAPS_WT}/s_${SAMPLE}/h1.fa
        if [ ! -s ${SAMPLE_FA} ]; then
            echo "ERROR: WT FASTA missing: ${SAMPLE_FA}" >&2; exit 1
        fi
        OUT_DIR=${HAPS_SV}/s_${SAMPLE}_sv_del_${SIZE}
        if [ -s ${OUT_DIR}/h1.fa ] && [ -s ${OUT_DIR}/h1.fa.fai ]; then
            echo "[$(date)] Reusing ${OUT_DIR}/h1.fa"
            continue
        fi
        rm -rf ${OUT_DIR}; mkdir -p ${OUT_DIR}
        echo "[$(date)] HACk DEL ${SIZE} (len=${LEN}) on ${SAMPLE}"
        VISOR HACk -g ${SAMPLE_FA} -b ${BED} -o ${OUT_DIR}
        if [ ! -s ${OUT_DIR}/h1.fa ]; then
            echo "ERROR: HACk produced no output for ${SAMPLE}" >&2; exit 1
        fi
    done
done

echo "[$(date)] DONE. SV haplotypes:"
for SIZE in "${!DEL_SIZES[@]}"; do
    for SAMPLE in "${SV_SAMPLES[@]}"; do
        ls -lh "${HAPS_SV}/s_${SAMPLE}_sv_del_${SIZE}/h1.fa" 2>/dev/null || true
    done
done
