#!/bin/bash
#SBATCH --job-name=apply_vcf
#SBATCH --output=logs/00_apply_vcf_%j.out
#SBATCH --error=logs/00_apply_vcf_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# =============================================================================
# 00_apply_vcf.sh
# Purpose: For each of the first N_SAMPLES ecotypes in the GrENET VCF,
#          apply their SNPs to Chr1.fa using bcftools consensus, producing
#          a per-sample FASTA that will serve as the "reference" for HACk.
#
# Requires: 00_prep_vcf.sh must have been run first (produces RENAMED_VCF
#           with chromosomes already named "Chr1").
#
# Key design decisions:
#  - Uses RENAMED_VCF (chrom "Chr1") → no inline renaming needed
#  - -H 1 (first allele): inbred lines are mostly homozygous
#  - SNPs only in VCF → SV BED coordinates remain valid in consensus FASTAs
#  - Output: ${HAPS_VAR}/s_<sample>/h1.fa  (one FASTA per ecotype)
# =============================================================================

set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONFIG_FILE=${1:-"$(dirname "$0")/config_sv_var_deletions.sh"}
source "${CONFIG_FILE}"

mkdir -p logs "${HAPS_VAR}"

# ---------------------------------------------------------------------------
# Fast early-exit: if all N_SAMPLES FASTAs already exist, nothing to do.
# This check runs BEFORE conda activation to avoid unnecessary overhead.
# ---------------------------------------------------------------------------
N_EXISTING=$(find "${HAPS_VAR}" -mindepth 2 -maxdepth 2 -name "h1.fa" 2>/dev/null | wc -l) || N_EXISTING=0
if [[ "${N_EXISTING}" -ge "${N_SAMPLES}" ]]; then
  echo "[$(date)] Skipping 00_apply_vcf — all ${N_SAMPLES} haplotype FASTAs already exist in ${HAPS_VAR}"
  echo "[$(date)] To force regeneration, remove the directory: rm -rf ${HAPS_VAR}"
  exit 0
fi
echo "[$(date)] Found ${N_EXISTING}/${N_SAMPLES} haplotypes — building missing ones."

source "$(mamba info --base)/etc/profile.d/conda.sh" && conda activate pang

# ---------------------------------------------------------------------------
# Guard: renamed VCF must exist
# ---------------------------------------------------------------------------
if [[ ! -s "${RENAMED_VCF}" || ! -s "${RENAMED_VCF}.tbi" ]]; then
  echo "ERROR: ${RENAMED_VCF} not found or not indexed." >&2
  echo "       Run 00_prep_vcf.sh first." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# Collect sample list (first N_SAMPLES from VCF header)
# ---------------------------------------------------------------------------
mapfile -t SAMPLES < <(bcftools query -l "${RENAMED_VCF}" | head -n "${N_SAMPLES}")
echo "[$(date)] Using ${#SAMPLES[@]} samples: ${SAMPLES[*]}"

# ---------------------------------------------------------------------------
# Per-sample bcftools consensus → Chr1 FASTA with sample SNPs applied
# ---------------------------------------------------------------------------
for SAMPLE in "${SAMPLES[@]}"; do
  OUT_DIR="${HAPS_VAR}/s_${SAMPLE}"
  OUT_FA="${OUT_DIR}/h1.fa"

  if [[ -s "${OUT_FA}" ]]; then
    echo "[$(date)] Skipping ${SAMPLE} — ${OUT_FA} already exists"
    continue
  fi

  mkdir -p "${OUT_DIR}"
  echo "[$(date)] Building consensus FASTA for sample ${SAMPLE}"

  # bcftools consensus prints "Applied N variants" to stderr — redirect to stdout
  # so the .err file stays empty on clean runs
  bcftools consensus \
      --fasta-ref "${REF}" \
      --samples  "${SAMPLE}" \
      --haplotype 1 \
      --output   "${OUT_FA}" \
      "${RENAMED_VCF}" 2>&1

  # Sanity check: exactly 1 sequence
  N_SEQS=$(grep -c "^>" "${OUT_FA}" || true)
  if [[ "${N_SEQS}" -ne 1 ]]; then
    echo "ERROR: ${OUT_FA} contains ${N_SEQS} sequences (expected 1)" >&2
    rm -f "${OUT_FA}"
    exit 1
  fi
  echo "[$(date)]   Done: ${OUT_FA}  ($(du -sh "${OUT_FA}" | cut -f1))"
done

echo "[$(date)] All ${#SAMPLES[@]} consensus FASTAs ready in ${HAPS_VAR}"
ls -lh "${HAPS_VAR}"
