#!/bin/bash
#SBATCH --job-name=visor_hack_var
#SBATCH --output=logs/02_run_hack_var_%j.out
#SBATCH --error=logs/02_run_hack_var_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# =============================================================================
# 02_run_hack_var.sh
# Purpose: Run VISOR HACk on the first N_SV = round(N_SAMPLES × SV_FREQ)
#          per-sample consensus FASTAs produced by 00_apply_vcf.sh.
#
# Each SV-bearing sample gets its own VISOR HACk output directory:
#   ${HAPS_VAR}/s_<sample>_sv/   ← contains h1.fa (consensus + deletion)
#
# Non-SV samples are left untouched in ${HAPS_VAR}/s_<sample>/
# (already created by 00_apply_vcf.sh).
#
# Key: because the GrENET VCF contains SNPs only (no indels), the deletion
#      BED coordinates are identical in the per-sample and original FASTAs.
# =============================================================================

set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
source "$(mamba info --base)/etc/profile.d/conda.sh" && conda activate pang

CONFIG_FILE=${1:-"$(dirname "$0")/config_sv_var_deletions.sh"}
source "${CONFIG_FILE}"

mkdir -p logs "${HAPS_VAR}"

# ---------------------------------------------------------------------------
# Helpers (same as 02_run_hack.sh)
# ---------------------------------------------------------------------------
validate_hap() {
  local fa="$1"
  local fai="${fa}.fai"
  if [[ ! -s "$fa" ]]; then
    # Expected on first run — informational, not an error
    echo "[validate] not found yet (will run HACk): $fa"
    return 1
  fi
  local n_seqs
  n_seqs=$(wc -l < "$fai" 2>/dev/null || echo 0)
  if [[ "$n_seqs" -ne 1 ]]; then
    # Genuine error — corrupted output from VISOR
    echo "[validate] ERROR: CORRUPTED: $fa has ${n_seqs} sequences in .fai (expected 1) — deleting" >&2
    rm -f "$fa" "$fai"
    return 1
  fi
  return 0
}

# ---------------------------------------------------------------------------
# Determine which samples exist (produced by 00_apply_vcf.sh) and how many
# carry the SV.
# ---------------------------------------------------------------------------
mapfile -t SAMPLES < <(bcftools query -l "${RENAMED_VCF}" | head -n "${N_SAMPLES}")

# N_SV = round(N_SAMPLES * SV_FREQ)
N_SV=$(python3 -c "import math; print(round(${N_SAMPLES} * ${SV_FREQ}))")
echo "[$(date)] N_SAMPLES=${N_SAMPLES}  SV_FREQ=${SV_FREQ}  N_SV=${N_SV}"

if [[ "${N_SV}" -lt 1 ]]; then
  echo "ERROR: N_SV=${N_SV} — increase N_SAMPLES or SV_FREQ" >&2
  exit 1
fi
if [[ "${N_SV}" -ge "${N_SAMPLES}" ]]; then
  echo "ERROR: N_SV=${N_SV} equals N_SAMPLES — all clones would carry the SV (freq=1)" >&2
  exit 1
fi

SV_SAMPLES=("${SAMPLES[@]:0:${N_SV}}")
echo "[$(date)] SV-bearing samples (${N_SV}): ${SV_SAMPLES[*]}"

# ---------------------------------------------------------------------------
# Run HACk for each SV size on each SV-bearing sample
# ---------------------------------------------------------------------------
case "${SV_TYPE}" in
  "DEL")
    for SIZE in "${!DEL_SIZES[@]}"; do
      LEN=${DEL_SIZES[$SIZE]}
      BED="${BEDS}/hack_del_${SIZE}.bed"

      if [[ ! -s "${BED}" ]]; then
        echo "ERROR: BED file not found: ${BED} — run 01_make_beds.sh first" >&2
        exit 1
      fi

      for SAMPLE in "${SV_SAMPLES[@]}"; do
        SAMPLE_FA="${HAPS_VAR}/s_${SAMPLE}/h1.fa"
        if [[ ! -s "${SAMPLE_FA}" ]]; then
          echo "ERROR: consensus FASTA missing: ${SAMPLE_FA} — run 00_apply_vcf.sh first" >&2
          exit 1
        fi

        OUT_DIR="${HAPS_VAR}/s_${SAMPLE}_sv_del_${SIZE}"

        if validate_hap "${OUT_DIR}/h1.fa"; then
          echo "[$(date)] Reusing existing SV haplotype: ${OUT_DIR}/h1.fa"
          continue
        fi

        rm -rf "${OUT_DIR}"
        mkdir -p "${OUT_DIR}"

        echo "[$(date)] HACk DEL ${SIZE} (len=${LEN}) on sample ${SAMPLE}"
        VISOR HACk -g "${SAMPLE_FA}" -b "${BED}" -o "${OUT_DIR}"

        validate_hap "${OUT_DIR}/h1.fa" \
          || { echo "ERROR: HACk produced corrupt FASTA for ${SAMPLE} DEL ${SIZE}" >&2; exit 1; }

        echo "[$(date)]   Done: ${OUT_DIR}/h1.fa"
      done
    done
    ;;
  *)
    echo "02_run_hack_var.sh: unsupported SV_TYPE=${SV_TYPE}" >&2
    exit 1
    ;;
esac

echo "[$(date)] All SV haplotypes created. Summary:"
for SIZE in "${!DEL_SIZES[@]}"; do
  for SAMPLE in "${SV_SAMPLES[@]}"; do
    ls -lh "${HAPS_VAR}/s_${SAMPLE}_sv_del_${SIZE}/h1.fa" 2>/dev/null || true
  done
done
