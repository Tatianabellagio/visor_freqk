#!/bin/bash
#SBATCH --job-name=prep_vcf
#SBATCH --output=logs/00_prep_vcf_%j.out
#SBATCH --error=logs/00_prep_vcf_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# =============================================================================
# 00_prep_vcf.sh  — run ONCE before anything else
# Purpose: Rename chromosome "1" → "Chr1" in the GrENET VCF so that all
#          downstream scripts (bcftools consensus, etc.) work directly against
#          the Chr1.fa reference without any inline renaming.
#
# Input:  ${VCF_FILE}          (original VCF, chrom = "1")
# Output: ${RENAMED_VCF}       (new VCF, chrom = "Chr1", bgzipped + tabix)
#         ${CHROM_MAP}         (the two-column mapping file used here)
#
# This script is idempotent: if RENAMED_VCF already exists and is indexed,
# it exits immediately.
# =============================================================================

set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
source "$(mamba info --base)/etc/profile.d/conda.sh" && conda activate pang

CONFIG_FILE=${1:-"$(dirname "$0")/config_sv_var_deletions.sh"}
source "${CONFIG_FILE}"

mkdir -p logs "$(dirname "${RENAMED_VCF}")"

# ---------------------------------------------------------------------------
# Fast exit if already done
# ---------------------------------------------------------------------------
if [[ -s "${RENAMED_VCF}" && -s "${RENAMED_VCF}.tbi" ]]; then
  echo "[$(date)] ${RENAMED_VCF} already exists and is indexed — nothing to do."
  bcftools index -s "${RENAMED_VCF}"
  exit 0
fi

# ---------------------------------------------------------------------------
# Write chrom map: "1" → "Chr1"  (tab-separated, one pair per line)
# ---------------------------------------------------------------------------
echo "[$(date)] Writing chrom map: ${CHROM_MAP}"
mkdir -p "$(dirname "${CHROM_MAP}")"
printf '1\tChr1\n' > "${CHROM_MAP}"

# ---------------------------------------------------------------------------
# Rename + compress + index
# Only chromosome 1 is needed downstream; use --regions 1 to skip 2-5
# and keep the output small.
# ---------------------------------------------------------------------------
echo "[$(date)] Renaming chromosomes in VCF (this may take a few minutes)..."
bcftools annotate \
    --rename-chrs "${CHROM_MAP}" \
    --regions 1 \
    --output-type z \
    --threads 4 \
    --output "${RENAMED_VCF}" \
    "${VCF_FILE}"

echo "[$(date)] Indexing ${RENAMED_VCF}..."
bcftools index --tbi --threads 4 "${RENAMED_VCF}"

echo "[$(date)] Done. Chromosome summary of renamed VCF:"
bcftools index -s "${RENAMED_VCF}"

echo "[$(date)] Sample count: $(bcftools query -l "${RENAMED_VCF}" | wc -l)"
echo "[$(date)] Renamed VCF ready: ${RENAMED_VCF}"
