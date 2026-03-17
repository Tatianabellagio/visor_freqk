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

CONFIG_FILE=${1:-"$(dirname "$0")/config_sv_var_deletions.sh"}
source "${CONFIG_FILE}"

mkdir -p logs "${HAPS_VAR}"

# ---------------------------------------------------------------------------
# Fast early-exit: if all SV haplotypes (N_SV × N_SIZES) already exist, skip.
# Runs before conda activation to avoid unnecessary overhead.
# ---------------------------------------------------------------------------
_N_SV=$(python3 -c "print(round(${N_SAMPLES} * ${SV_FREQ}))")
_N_SIZES=${#DEL_SIZES[@]}
_N_HAP_EXPECTED=$((_N_SV * _N_SIZES))
_N_HAP_EXISTING=$(find "${HAPS_VAR}" -mindepth 2 -maxdepth 2 -name "h1.fa" \
                    -path "*_sv_del_*" -size +0c 2>/dev/null | wc -l) || _N_HAP_EXISTING=0
if [[ "${_N_HAP_EXPECTED}" -gt 0 && "${_N_HAP_EXISTING}" -ge "${_N_HAP_EXPECTED}" ]]; then
  echo "[$(date)] Skipping 02_run_hack_var — all ${_N_HAP_EXPECTED} SV haplotypes already exist in ${HAPS_VAR}"
  echo "[$(date)] To force regeneration: rm -rf ${HAPS_VAR}/*_sv_del_*"
  exit 0
fi
echo "[$(date)] Found ${_N_HAP_EXISTING}/${_N_HAP_EXPECTED} SV haplotypes — building missing ones."

source "$(mamba info --base)/etc/profile.d/conda.sh" && conda activate pang

# ---------------------------------------------------------------------------
# Helpers (same as 02_run_hack.sh)
# ---------------------------------------------------------------------------
validate_hap() {
  local fa="$1"
  local fai="${fa}.fai"
  if [[ ! -s "$fa" ]]; then
    echo "[validate] not found yet (will run HACk): $fa"
    return 1
  fi
  local n_seqs
  n_seqs=$(wc -l < "$fai" 2>/dev/null || echo 0)
  if [[ "$n_seqs" -ne 1 ]]; then
    echo "[validate] ERROR: CORRUPTED: $fa has ${n_seqs} sequences in .fai (expected 1) — deleting" >&2
    rm -f "$fa" "$fai" 2>/dev/null || true   # NFS: ignore busy-file errors
    return 1
  fi
  return 0
}

# ---------------------------------------------------------------------------
# Per-sample NFS-safe locking.
#   acquire_lock DIR   — blocks until lock is owned or times out (exit 1)
#   release_lock DIR   — removes the lock; safe to call even if not held
#
# mkdir is atomic on NFS, so two concurrent callers cannot both succeed.
# If the lock dir is older than MAX_LOCK_AGE_SECS we treat it as stale
# (the job that created it was likely killed by Slurm) and steal it.
# ---------------------------------------------------------------------------
MAX_LOCK_AGE_SECS=7200   # 2 hours — longer than the job time limit

acquire_lock() {
  local out_dir="$1"
  local lock="${out_dir}.lock"
  local waited=0
  local poll=15          # seconds between retries

  while true; do
    if mkdir "${lock}" 2>/dev/null; then
      # We own the lock — register cleanup so Slurm kills don't leave stale locks
      _CURRENT_LOCK="${lock}"
      return 0
    fi

    # Lock exists — check if it's stale
    if [[ -d "${lock}" ]]; then
      local age
      age=$(( $(date +%s) - $(stat -c %Y "${lock}" 2>/dev/null || echo $(date +%s)) ))
      if [[ "${age}" -gt "${MAX_LOCK_AGE_SECS}" ]]; then
        echo "[lock] WARNING: stale lock ${lock} (age=${age}s) — stealing it" >&2
        rmdir "${lock}" 2>/dev/null || true
        continue   # retry mkdir
      fi
    fi

    # Fresh lock — wait
    waited=$(( waited + poll ))
    if [[ "${waited}" -ge 3600 ]]; then
      echo "ERROR: timed out waiting for lock on ${out_dir} after ${waited}s" >&2
      exit 1
    fi
    echo "[lock] waiting for ${lock} (${waited}s elapsed)..."
    sleep "${poll}"

    # If the h1.fa appeared while we were waiting the other job succeeded — done
    if [[ -s "${out_dir}/h1.fa" ]]; then
      _CURRENT_LOCK=""
      return 2   # sentinel: already built by concurrent job
    fi
  done
}

release_lock() {
  local lock="$1"
  rmdir "${lock}" 2>/dev/null || true
  _CURRENT_LOCK=""
}

_CURRENT_LOCK=""
trap 'release_lock "${_CURRENT_LOCK}" 2>/dev/null || true' EXIT

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

        # Fast path — already built and valid
        if validate_hap "${OUT_DIR}/h1.fa"; then
          echo "[$(date)] Reusing existing SV haplotype: ${OUT_DIR}/h1.fa"
          continue
        fi

        # Acquire per-sample lock (NFS-safe atomic mkdir).
        # Returns 0 = lock owned, 2 = concurrent job built it while we waited.
        acquire_lock "${OUT_DIR}"
        lock_rc=$?

        if [[ "${lock_rc}" -eq 2 ]]; then
          # Another job finished while we were waiting — validate and move on
          validate_hap "${OUT_DIR}/h1.fa" \
            || { echo "ERROR: concurrent build of ${SAMPLE} DEL ${SIZE} left corrupt FASTA" >&2; exit 1; }
          echo "[$(date)] Built by concurrent job: ${OUT_DIR}/h1.fa"
          continue
        fi

        # We own the lock — clean up any corrupt partial output without
        # rm -rf (NFS ghost files (.nfs*) make rm -rf fail on busy dirs).
        rm -f "${OUT_DIR}/h1.fa" "${OUT_DIR}/h1.fa.fai" 2>/dev/null || true
        mkdir -p "${OUT_DIR}"

        echo "[$(date)] HACk DEL ${SIZE} (len=${LEN}) on sample ${SAMPLE}"
        VISOR HACk -g "${SAMPLE_FA}" -b "${BED}" -o "${OUT_DIR}"

        if validate_hap "${OUT_DIR}/h1.fa"; then
          echo "[$(date)]   Done: ${OUT_DIR}/h1.fa"
          release_lock "${OUT_DIR}.lock"
        else
          release_lock "${OUT_DIR}.lock"
          echo "ERROR: HACk produced corrupt FASTA for ${SAMPLE} DEL ${SIZE}" >&2
          exit 1
        fi
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
