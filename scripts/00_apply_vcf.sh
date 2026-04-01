#!/bin/bash
#SBATCH --job-name=apply_vcf
#SBATCH --output=logs/00_apply_vcf_%j.out
#SBATCH --error=logs/00_apply_vcf_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# =============================================================================
# 00_apply_vcf.sh
# Purpose: For each of the first N_SAMPLES ecotypes in the GrENET VCF,
#          apply their SNPs to Chr1.fa using bcftools consensus, producing
#          a per-sample FASTA that serves as the "reference" for HACk.
#
# Output goes to HAPS_WT (shared across all positions, keyed by n_samples):
#   ${HAPS_WT}/s_<sample>/h1.fa
#
# Because HAPS_WT is shared, multiple pipelines for different positions may
# run 00_apply_vcf simultaneously.  A flock-based check-lock-check pattern
# ensures only one process builds the FASTAs; others wait then skip.
# Reading finished FASTAs from multiple pipelines in parallel is safe.
# =============================================================================

set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONFIG_FILE=${1:-"$(dirname "$0")/config_sv_var_deletions.sh"}
source "${CONFIG_FILE}"

mkdir -p logs "${HAPS_WT}"

LOCK_FILE="${HAPS_WT}/.apply_vcf.lock"

_n_existing() {
  find "${HAPS_WT}" -mindepth 2 -maxdepth 2 -name "h1.fa" 2>/dev/null | wc -l || echo 0
}

# ---------------------------------------------------------------------------
# Fast path (no lock): all FASTAs already exist
# ---------------------------------------------------------------------------
N_EXISTING=$(_n_existing)
if [[ "${N_EXISTING}" -ge "${N_SAMPLES}" ]]; then
  echo "[$(date)] Skipping 00_apply_vcf — all ${N_SAMPLES} WT FASTAs already exist in ${HAPS_WT}"
  exit 0
fi

echo "[$(date)] Found ${N_EXISTING}/${N_SAMPLES} WT FASTAs — acquiring build lock …"

# ---------------------------------------------------------------------------
# Slow path: exclusive flock → re-check → build.
# A concurrent job may have finished while we waited for the lock.
# ---------------------------------------------------------------------------
(
  flock -x 200

  N_EXISTING=$(_n_existing)
  if [[ "${N_EXISTING}" -ge "${N_SAMPLES}" ]]; then
    echo "[$(date)] WT FASTAs built by concurrent job — nothing to do."
    exit 0
  fi
  echo "[$(date)] Lock acquired. Building ${N_SAMPLES} WT FASTAs in ${HAPS_WT}"

  source "$(mamba info --base)/etc/profile.d/conda.sh" && conda activate pang

  if [[ ! -s "${RENAMED_VCF}" || ! -s "${RENAMED_VCF}.tbi" ]]; then
    echo "ERROR: ${RENAMED_VCF} not found or not indexed — run 00_prep_vcf.sh first." >&2
    exit 1
  fi

  mapfile -t SAMPLES < <(bcftools query -l "${RENAMED_VCF}" | head -n "${N_SAMPLES}")
  echo "[$(date)] Using ${#SAMPLES[@]} samples: ${SAMPLES[*]}"

  for SAMPLE in "${SAMPLES[@]}"; do
    OUT_DIR="${HAPS_WT}/s_${SAMPLE}"
    OUT_FA="${OUT_DIR}/h1.fa"

    if [[ -s "${OUT_FA}" ]]; then
      echo "[$(date)] Skipping ${SAMPLE} — already exists"
      continue
    fi

    mkdir -p "${OUT_DIR}"
    echo "[$(date)] Building consensus FASTA for ${SAMPLE}"

    bcftools consensus \
        --fasta-ref "${REF}" \
        --samples   "${SAMPLE}" \
        --haplotype 1 \
        --output    "${OUT_FA}" \
        "${RENAMED_VCF}" 2>&1

    N_SEQS=$(grep -c "^>" "${OUT_FA}" || true)
    if [[ "${N_SEQS}" -ne 1 ]]; then
      echo "ERROR: ${OUT_FA} has ${N_SEQS} sequences (expected 1)" >&2
      rm -f "${OUT_FA}"
      exit 1
    fi
    echo "[$(date)]   Done: ${OUT_FA}  ($(du -sh "${OUT_FA}" | cut -f1))"
  done

  echo "[$(date)] All ${#SAMPLES[@]} WT FASTAs ready in ${HAPS_WT}"

) 200>"${LOCK_FILE}"
