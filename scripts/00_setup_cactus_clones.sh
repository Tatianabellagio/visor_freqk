#!/bin/bash
#SBATCH --job-name=cactus_setup
#SBATCH --output=logs/00_setup_cactus_clones_%j.out
#SBATCH --error=logs/00_setup_cactus_clones_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G

# =============================================================================
# 00_setup_cactus_clones.sh
# Replaces 00_apply_vcf for the cactus_em pipeline.
#
# Picks N_SAMPLES cactus founders and creates the HAPS_WT/s_<id>/h1.fa
# structure as symlinks to the cactus FASTAs. Also writes the SHORtS
# region BED restricting to CHROM (Chr1 by default).
# =============================================================================

set -euo pipefail
CONFIG_FILE=${1:-"$(dirname "$0")/config_cactus_em.sh"}
source "${CONFIG_FILE}"

mkdir -p logs "${HAPS_WT}" "${BEDS}"

# List cactus founders alphabetically, pick first N (deterministic)
mapfile -t FOUNDERS < <(ls ${CACTUS_DIR}/*.chr.fa 2>/dev/null | \
    xargs -n1 basename | sed 's/\.chr\.fa$//' | grep -v "^TAIR10$" | sort | head -n ${N_SAMPLES})

if [ ${#FOUNDERS[@]} -lt ${N_SAMPLES} ]; then
    echo "ERROR: only ${#FOUNDERS[@]} cactus founders available (need ${N_SAMPLES})" >&2
    exit 1
fi
echo "[$(date)] Using ${#FOUNDERS[@]} cactus founders: ${FOUNDERS[*]}"

# Link each founder's FASTA into HAPS_WT/s_<id>/h1.fa
for FID in "${FOUNDERS[@]}"; do
    OUT_DIR=${HAPS_WT}/s_${FID}
    mkdir -p ${OUT_DIR}
    ln -sfn ${CACTUS_DIR}/${FID}.chr.fa     ${OUT_DIR}/h1.fa
    ln -sfn ${CACTUS_DIR}/${FID}.chr.fa.fai ${OUT_DIR}/h1.fa.fai
done
echo "[$(date)] Linked ${#FOUNDERS[@]} WT haplotypes in ${HAPS_WT}"

# Build SHORtS region BED (chrom-wide, format: chrom\tstart\tend\tcap_pct)
SAMTOOLS=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/samtools
${SAMTOOLS} faidx ${REF} 2>/dev/null || true
CHR_LEN=$(awk -v c=${CHROM} '$1==c {print $2}' ${REF}.fai)
if [ -z "${CHR_LEN}" ]; then
    echo "ERROR: chrom ${CHROM} not found in ${REF}.fai" >&2; exit 1
fi
echo -e "${CHROM}\t0\t${CHR_LEN}\t100.0" > ${REGION_BED_PATH}
echo "[$(date)] SHORtS region BED:"
cat ${REGION_BED_PATH}

# Save sample list for downstream
printf "%s\n" "${FOUNDERS[@]}" > ${HAPS_WT}/../samples.txt
echo "[$(date)] Sample list: ${HAPS_WT}/../samples.txt"
