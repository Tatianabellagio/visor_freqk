#!/bin/bash -l
#SBATCH --job-name=visor_shorts_var
#SBATCH --output=logs/03_run_shorts_var_%j.out
#SBATCH --error=logs/03_run_shorts_var_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00

# =============================================================================
# 03_run_shorts_var.sh
# Purpose: Simulate pool-seq reads from N ecotype haplotypes where
#          N_SV = round(N × SV_FREQ) carry a deletion and the rest are WT.
#
#   Clone list passed to SHORtS:
#     s_<sample>_sv_del_<SIZE>/  × N_SV   (deletion haplotypes)
#     s_<sample>/                × N-N_SV  (WT haplotypes with sample SNPs)
#
#   Each clone contributes equal coverage: clonefraction = 100/N each.
#   The last fraction absorbs any floating-point remainder so fracs sum
#   to exactly 100 (e.g. N=3 → 33.3333 33.3333 33.3334).
#
# Note: SHORtS accepts any number of -s dirs; there is no hardcoded limit.
# =============================================================================

set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONFIG_FILE=${1:-"$(dirname "$0")/config_sv_var_deletions.sh"}
source "${CONFIG_FILE}"

mkdir -p logs "${BEDS}"

# ---------------------------------------------------------------------------
# Fast early-exit: if all read sets already exist, skip conda + SHORtS.
# Runs before conda activation to avoid unnecessary overhead.
# ---------------------------------------------------------------------------
_FL=$(python3 -c "print(round(${SV_FREQ}*100))")
_EL=$(python3 -c "e='${ERROR_RATE}'; print('0' if float(e)==0 else (e.split('.')[-1] if '.' in e else e))")
_RTAG="n${N_SAMPLES}_f${_FL}_err${_EL}"
_N_READS_DONE=0; _N_READS_TOTAL=0
for _SZ in "${!DEL_SIZES[@]}"; do
  _N_READS_TOTAL=$((_N_READS_TOTAL + 1))
  _OUTD="${READS_VAR}/cov${COVERAGE}/var_del_${_SZ}_${_RTAG}"
  [[ -s "${_OUTD}/r1.fq" && -s "${_OUTD}/r2.fq" ]] && _N_READS_DONE=$((_N_READS_DONE + 1))
done
if [[ "${_N_READS_TOTAL}" -gt 0 && "${_N_READS_DONE}" -ge "${_N_READS_TOTAL}" ]]; then
  echo "[$(date)] Skipping 03_run_shorts_var — all ${_N_READS_TOTAL} read sets already exist"
  echo "[$(date)] RUN_TAG=${_RTAG}  READS_VAR=${READS_VAR}/cov${COVERAGE}/"
  echo "[$(date)] To force regeneration, remove the read directories."
  exit 0
fi
echo "[$(date)] Found ${_N_READS_DONE}/${_N_READS_TOTAL} read sets — building missing ones."

source "$(mamba info --base)/etc/profile.d/conda.sh" && conda activate pang

# ---------------------------------------------------------------------------
# Sample list and SV assignments (must match 02_run_hack_var.sh)
# ---------------------------------------------------------------------------
mapfile -t SAMPLES < <(bcftools query -l "${RENAMED_VCF}" | head -n "${N_SAMPLES}")
N_SV=$(python3 -c "print(round(${N_SAMPLES} * ${SV_FREQ}))")
N_WT=$(( N_SAMPLES - N_SV ))
echo "[$(date)] N_SAMPLES=${N_SAMPLES}  N_SV=${N_SV}  N_WT=${N_WT}  SV_FREQ=${SV_FREQ}"

# ---------------------------------------------------------------------------
# Helper: compute equal clone fractions for N clones that sum to exactly 100.
# Fractions are derived from the actual CLONE_DIRS count (not hardcoded from
# N_SAMPLES), so missing haplotypes are handled gracefully.
# Avoids nested double-quotes by storing the comma-separated list in a variable.
# ---------------------------------------------------------------------------
make_fracs() {
  local n="$1"
  python3 - "${n}" << 'PY'
import sys
n = int(sys.argv[1])
base = round(100.0 / n, 4)
fracs = [base] * n
# Last fraction = 100 - sum(first n-1), using the same iterative sum
# VISOR will do.  This guarantees sum(fracs) == 100.0 exactly in float
# arithmetic — no decimal-rounding artefacts.
fracs[-1] = 100.0 - sum(fracs[:-1])
for f in fracs:
    print(f)
PY
}

# ---------------------------------------------------------------------------
# Region BED (full Chr1)
# ---------------------------------------------------------------------------
CHR1_LEN=$(awk -v c="${CHROM}" '$1==c {print $2}' "${REF}.fai")
REGION_BED="${BEDS}/shorts_region_${CHROM}.bed"
echo -e "${CHROM}\t1\t${CHR1_LEN}\t100.0\t100.0" > "${REGION_BED}"
echo "[$(date)] Region BED: ${CHROM} 1-${CHR1_LEN}"

# ---------------------------------------------------------------------------
# Run for each SV size
# ---------------------------------------------------------------------------
case "${SV_TYPE}" in
  "DEL")
    for SIZE in "${!DEL_SIZES[@]}"; do

      # Build ordered clone-dir list: SV first, then WT
      # (order must match the FRACS array built above)
      CLONE_DIRS=()

      for SAMPLE in "${SAMPLES[@]:0:${N_SV}}"; do
        SV_DIR="${HAPS_SV}/s_${SAMPLE}_sv_del_${SIZE}"
        if [[ ! -s "${SV_DIR}/h1.fa" ]]; then
          echo "ERROR: SV haplotype missing: ${SV_DIR}/h1.fa — run 02_run_hack_var.sh first" >&2
          exit 1
        fi
        CLONE_DIRS+=("${SV_DIR}")
      done

      for SAMPLE in "${SAMPLES[@]:${N_SV}}"; do
        WT_DIR="${HAPS_WT}/s_${SAMPLE}"
        if [[ ! -s "${WT_DIR}/h1.fa" ]]; then
          echo "ERROR: WT haplotype missing: ${WT_DIR}/h1.fa — run 00_apply_vcf.sh first" >&2
          exit 1
        fi
        CLONE_DIRS+=("${WT_DIR}")
      done

      # Compute fractions from the actual number of clones assembled above.
      # This handles any mismatch gracefully (e.g. a missing haplotype skipped).
      N_CLONES=${#CLONE_DIRS[@]}
      mapfile -t FRACS < <(make_fracs "${N_CLONES}")
      printf -v FRAC_CSV '%s,' "${FRACS[@]}"; FRAC_CSV="${FRAC_CSV%,}"  # join with commas
      FRAC_SUM=$(python3 -c "print(round(sum([${FRAC_CSV}]),4))")
      echo "[$(date)] Clones: ${N_CLONES}  fractions: ${FRACS[*]}  sum=${FRAC_SUM}"
      echo "[$(date)]   SV clones (${N_SV}): ${CLONE_DIRS[*]:0:${N_SV}}"
      echo "[$(date)]   WT clones (${N_WT}): ${CLONE_DIRS[*]:${N_SV}}"

      # Encode labels for output directory
      FREQ_LABEL=$(python3 -c "print(round(${SV_FREQ}*100))")
      ERR_LABEL=$(python3 -c "
e='${ERROR_RATE}'
if float(e)==0: print('0')
else: print(e.split('.')[-1] if '.' in e else e)
")
      RUN_TAG="n${N_SAMPLES}_f${FREQ_LABEL}_err${ERR_LABEL}"
      OUT="${READS_VAR}/cov${COVERAGE}/var_del_${SIZE}_${RUN_TAG}"

      if [[ -s "${OUT}/r1.fq" && -s "${OUT}/r2.fq" ]]; then
        echo "[$(date)] Skipping SHORtS DEL ${SIZE} — reads already exist in ${OUT}"
        continue
      fi

      # Clean up any partial state from a previous failed run.
      # VISOR SHORtS uses os.rmdir() to clean up its clone*/h*/ temp dirs;
      # if a previous run was killed those dirs are non-empty and cause an
      # immediate OSError on the next attempt.  Wiping the whole output dir
      # (safe here because r1.fq/r2.fq don't exist yet) avoids this.
      rm -rf "${OUT}"
      mkdir -p "${OUT}"

      echo "[$(date)] Running VISOR SHORtS: ${N_CLONES} clones, DEL ${SIZE}, cov=${COVERAGE}x"

      # VISOR SHORtS has a bug: os.rmdir() on its own clone*/h*/ temp dirs fails
      # with "Directory not empty" even though VISOR created those files itself.
      # The reads ARE produced correctly before the cleanup step, so we use || true
      # to ignore the spurious non-zero exit, then validate r1.fq/r2.fq below.
      VISOR SHORtS \
          -g "${REF}" \
          -s "${CLONE_DIRS[@]}" \
          -b "${REGION_BED}" \
          -o "${OUT}" \
          --coverage "${COVERAGE}" \
          --clonefraction "${FRACS[@]}" \
          --error "${ERROR_RATE}" \
          --fastq \
          --threads 8 || true

      if [[ -s "${OUT}/r1.fq" && -s "${OUT}/r2.fq" ]]; then
        echo "[$(date)] Done. Output in ${OUT}"
        ls -lh "${OUT}"
      else
        echo "ERROR: SHORtS finished but r1.fq/r2.fq missing in ${OUT}" >&2
        exit 1
      fi
    done
    ;;
  *)
    echo "03_run_shorts_var.sh: unsupported SV_TYPE=${SV_TYPE}" >&2
    exit 1
    ;;
esac
