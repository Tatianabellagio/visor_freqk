#!/bin/bash
#SBATCH --job-name=visor_hack
#SBATCH --output=logs/02_run_hack_%j.out
#SBATCH --error=logs/02_run_hack_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
# =============================================================================
# 02_run_hack.sh
# Purpose: Run VISOR HACk to inject homozygous SVs into Chr1
#          (deletions or insertions, depending on config SV_TYPE)
#          Parameters come from a config file (default: config_sv_deletions.sh)
# Output (DEL): ${HAPS}/del_<size>/HAP1/ and HAP2/
# Output (INS): ${HAPS}/ins_<size>/HAP1/ and HAP2/
# =============================================================================
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
source "$(mamba info --base)/etc/profile.d/conda.sh" && conda activate pang

CONFIG_FILE=${1:-"$(dirname "$0")/config_sv_deletions.sh"}
source "${CONFIG_FILE}"

mkdir -p logs "${HAPS}"

# -----------------------------------------------------------------------------
# validate_hap <fa_path>
# Returns 0 if the haplotype FASTA is healthy (exactly 1 sequence in .fai).
# If corrupted (duplicate sequence headers, wrong .fai line count), prints a
# warning, deletes the bad fa + fai, and returns 1 so the caller re-runs HACk.
# Using .fai line count is fast (tiny file) and catches the VISOR HACk bug
# that occasionally writes Chr1 twice, producing a ~2x-sized FASTA.
# -----------------------------------------------------------------------------
validate_hap() {
  local fa="$1"
  local fai="${fa}.fai"
  if [[ ! -s "$fa" ]]; then
    echo "[validate] MISSING or empty: $fa" >&2
    return 1
  fi
  local n_seqs
  n_seqs=$(wc -l < "$fai" 2>/dev/null || echo 0)
  if [[ "$n_seqs" -ne 1 ]]; then
    echo "[validate] CORRUPTED: $fa has ${n_seqs} sequence(s) in .fai (expected 1) — deleting" >&2
    rm -f "$fa" "$fai"
    return 1
  fi
  return 0
}

# -----------------------------------------------------------------------------
# acquire_lock <lock_dir> [max_wait_seconds]
# NFS-safe mutex using atomic mkdir. Spins until the lock directory is created
# by this process, then registers a trap to release it on exit.
# Multiple concurrent 02_run_hack jobs share haplotype paths (same sv_type +
# pos_label), so without locking they race to write the same h1.fa on NFS,
# producing duplicate Chr1 headers and "Stale file handle" errors.
# -----------------------------------------------------------------------------
acquire_lock() {
  local lock="$1"
  local max_wait="${2:-300}"
  local waited=0
  while ! mkdir "$lock" 2>/dev/null; do
    echo "[lock] Waiting for $lock (held by another job)..."
    sleep $((RANDOM % 10 + 5))
    waited=$((waited + 10))
    if [[ $waited -ge $max_wait ]]; then
      echo "[lock] Timeout waiting for $lock — removing stale lock and retrying" >&2
      rm -rf "$lock"
    fi
  done
  # Release lock on any exit (normal, error, or signal)
  trap "rm -rf '${lock}'" EXIT INT TERM
}

release_lock() {
  local lock="$1"
  rm -rf "$lock"
  trap - EXIT INT TERM
}

# -----------------------------------------------------------------------------
# Create WT clone (reference with no variants) for SHORtS. SHORtS expects each
# sample dir to contain *.fa; 03_run_shorts uses WT_CLONE as the second clone.
# Without this, SHORtS would fail (or use wrong/missing path).
# -----------------------------------------------------------------------------
WT_CLONE_DIR="${WT_CLONE}"
mkdir -p "${WT_CLONE_DIR}"
REF_FA="${WT_CLONE_DIR}/$(basename "${REF}")"
if [[ ! -s "${REF_FA}" ]]; then
  echo "[$(date)] Creating WT clone at ${WT_CLONE_DIR} (copy of reference)"
  cp "${REF}" "${REF_FA}"
  cp "${REF}.fai" "${WT_CLONE_DIR}/$(basename "${REF}.fai")"
  echo "[$(date)] WT clone ready: ${REF_FA}"
else
  echo "[$(date)] WT clone already exists: ${REF_FA}"
fi

case "${SV_TYPE}" in
  "DEL")
    for SIZE in "${!DEL_SIZES[@]}"; do
        LEN=${DEL_SIZES[$SIZE]}
        BED=${BEDS}/hack_del_${SIZE}.bed
        HAP1_OUT=${HAPS}/del_${SIZE}/HAP1
        HAP2_OUT=${HAPS}/del_${SIZE}/HAP2

        LOCK="${HAPS}/del_${SIZE}.lock"
        acquire_lock "$LOCK"

        # If haplotypes already exist AND are valid (1 sequence each), skip recomputation
        if validate_hap "${HAP1_OUT}/h1.fa" && validate_hap "${HAP2_OUT}/h1.fa"; then
          echo "[$(date)] Reusing existing haplotypes for DEL ${SIZE} in ${HAP1_OUT}, ${HAP2_OUT}"
        else
          rm -rf "${HAP1_OUT}" "${HAP2_OUT}"
          mkdir -p "${HAP1_OUT}" "${HAP2_OUT}"

          echo "[$(date)] Running VISOR HACk (DEL) for size: ${SIZE} (len=${LEN})"

          VISOR HACk -g "${REF}" -b "${BED}" -o "${HAP1_OUT}"
          validate_hap "${HAP1_OUT}/h1.fa" || { echo "HACk produced corrupt HAP1 for DEL ${SIZE}" >&2; exit 1; }

          VISOR HACk -g "${REF}" -b "${BED}" -o "${HAP2_OUT}"
          validate_hap "${HAP2_OUT}/h1.fa" || { echo "HACk produced corrupt HAP2 for DEL ${SIZE}" >&2; exit 1; }

          echo "[$(date)] Done: del_${SIZE}"
        fi

        release_lock "$LOCK"
    done
    echo "[$(date)] All deletions processed. Haplotypes in ${HAPS}"
    ;;
  "INS")
    for SIZE in "${!INS_SIZES[@]}"; do
        LEN=${INS_SIZES[$SIZE]}
        BED=${BEDS}/hack_ins_${SIZE}.bed
        HAP1_OUT=${HAPS}/ins_${SIZE}/HAP1
        HAP2_OUT=${HAPS}/ins_${SIZE}/HAP2

        LOCK="${HAPS}/ins_${SIZE}.lock"
        acquire_lock "$LOCK"

        if validate_hap "${HAP1_OUT}/h1.fa" && validate_hap "${HAP2_OUT}/h1.fa"; then
          echo "[$(date)] Reusing existing haplotypes for INS ${SIZE} in ${HAP1_OUT}, ${HAP2_OUT}"
        else
          rm -rf "${HAP1_OUT}" "${HAP2_OUT}"
          mkdir -p "${HAP1_OUT}" "${HAP2_OUT}"

          echo "[$(date)] Running VISOR HACk (INS) for size: ${SIZE} (len=${LEN})"

          VISOR HACk -g "${REF}" -b "${BED}" -o "${HAP1_OUT}"
          validate_hap "${HAP1_OUT}/h1.fa" || { echo "HACk produced corrupt HAP1 for INS ${SIZE}" >&2; exit 1; }

          VISOR HACk -g "${REF}" -b "${BED}" -o "${HAP2_OUT}"
          validate_hap "${HAP2_OUT}/h1.fa" || { echo "HACk produced corrupt HAP2 for INS ${SIZE}" >&2; exit 1; }

          echo "[$(date)] Done: ins_${SIZE}"
        fi

        release_lock "$LOCK"
    done
    echo "[$(date)] All insertions processed. Haplotypes in ${HAPS}"
    ;;
  *)
    echo "02_run_hack.sh: unsupported SV_TYPE=${SV_TYPE} (expected DEL or INS)" >&2
    exit 1
    ;;
esac