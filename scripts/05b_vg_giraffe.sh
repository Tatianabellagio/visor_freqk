#!/bin/bash -l
#SBATCH --job-name=vg_giraffe
#SBATCH --output=logs/05b_vg_giraffe_%j.out
#SBATCH --error=logs/05b_vg_giraffe_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
# =============================================================================
# 05b_vg_giraffe.sh
# Purpose: Benchmark SV allele-frequency estimation with vg giraffe + vg call.
#          Reuses the same reads and truth VCF produced by 03_run_shorts_var.sh
#          and 04_make_vcf_var.sh.  No VISOR re-simulation.
#
#   graph index (cached per pos_label × SIZE):
#       data/vg_indexes/del/<pos>/<SIZE>/graph.{giraffe.gbz,dist,min[,zipcodes]}
#
#   per-run outputs (sibling of the freqk output TSV):
#       results/del/<pos>/var/cov<COV>_err<ERR>/<SIZE>/n<N>/f<F>/k<K>/
#           var_del_<SIZE>_<RUN_TAG>.vg_giraffe.aln.gam
#           var_del_<SIZE>_<RUN_TAG>.vg_giraffe.pack
#           var_del_<SIZE>_<RUN_TAG>.vg_giraffe.calls.vcf
#           var_del_<SIZE>_<RUN_TAG>.vg_giraffe.allele_frequencies.k<K>.tsv
#
# Output TSV schema matches freqk's: a single line "af_ref|af_alt".
#
# Usage:
#   sbatch scripts/05b_vg_giraffe.sh [config_file]
#   (default config: config_sv_var_deletions.sh)
# =============================================================================
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONFIG_FILE=${1:-"$(dirname "$0")/config_sv_var_deletions.sh"}
source "${CONFIG_FILE}"

# ---------------------------------------------------------------------------
# Labels (must mirror 03_run_shorts_var.sh and 05_freqk_var.sh exactly)
# ---------------------------------------------------------------------------
FREQ_LABEL=$(python3 -c "print(round(${SV_FREQ}*100))")
ERR_LABEL=$(python3 -c "
e='${ERROR_RATE}'
if float(e)==0: print('0')
else: print(e.split('.')[-1] if '.' in e else e)
")
COV_LABEL="cov${COVERAGE}"
ERR_LABEL_DIR="err${ERR_LABEL}"
FREQ_LABEL_DIR="f${FREQ_LABEL}"
N_LABEL="n${N_SAMPLES}"
K_LABEL="k${K}"
RUN_TAG="n${N_SAMPLES}_f${FREQ_LABEL}_err${ERR_LABEL}"

# ---------------------------------------------------------------------------
# Fast early-exit: skip if all per-size AF TSVs already exist
# ---------------------------------------------------------------------------
_N_DONE=0; _N_TOTAL=0
for _SZ in "${!DEL_SIZES[@]}"; do
  _N_TOTAL=$((_N_TOTAL + 1))
  _BDIR="${RESULTS}/${COV_LABEL}_${ERR_LABEL_DIR}/${_SZ}/${N_LABEL}/${FREQ_LABEL_DIR}/${K_LABEL}"
  _AF="${_BDIR}/var_del_${_SZ}_${RUN_TAG}.vg_giraffe.allele_frequencies.k${K}.tsv"
  [[ -s "${_AF}" ]] && _N_DONE=$((_N_DONE + 1))
done
if [[ "${_N_TOTAL}" -gt 0 && "${_N_DONE}" -ge "${_N_TOTAL}" ]]; then
  echo "[$(date)] Skipping 05b_vg_giraffe — all ${_N_TOTAL} AF files already exist"
  exit 0
fi
echo "[$(date)] Found ${_N_DONE}/${_N_TOTAL} AF files — running missing ones."

# ---------------------------------------------------------------------------
# Per-job scratch for vg's SDSL temp files.  vg autoindex / vg pack drop
# `text_*.sdsl`, `sa_*.sdsl`, `<pid>_*_text_rec*.sdsl`, etc. into cwd while
# building succinct data structures.  Without this, those land in the
# project root on NFS (slow + survive crashes as orphans).  cd into a
# per-job /tmp dir so they live on local disk and disappear with the job.
# All other paths used downstream are absolute (built from $WORK), so the
# cd is safe.
# ---------------------------------------------------------------------------
VG_SCRATCH="${SLURM_TMPDIR:-/tmp}/vg_${SLURM_JOB_ID:-$$}"
mkdir -p "${VG_SCRATCH}"
trap 'rm -rf "${VG_SCRATCH}"' EXIT
cd "${VG_SCRATCH}"
echo "[$(date)] vg scratch: ${VG_SCRATCH}"

# ---------------------------------------------------------------------------
# Env: pang has vg, samtools, bcftools (mosdepth not needed here)
# ---------------------------------------------------------------------------
source "$(mamba info --base)/etc/profile.d/conda.sh" && mamba activate pang

command -v vg       >/dev/null || { echo "ERROR: vg not found on PATH"       >&2; exit 1; }
command -v bcftools >/dev/null || { echo "ERROR: bcftools not found on PATH" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Index cache root (shared across cov/freq/err within the same pos × SIZE)
# ---------------------------------------------------------------------------
VG_INDEX_ROOT=${WORK}/data/vg_indexes/del/${POS_LABEL}
mkdir -p "${VG_INDEX_ROOT}"

# ---------------------------------------------------------------------------
# Helper: build vg giraffe index for one (pos_label, SIZE) if not cached.
# Uses a per-size lock so concurrent sv_freq jobs don't race.
# ---------------------------------------------------------------------------
ensure_index() {
  local size="$1" vcf="$2" outdir="$3"
  local gbz="${outdir}/graph.giraffe.gbz"
  local dist="${outdir}/graph.dist"
  # vg 1.70 autoindex emits graph.shortread.withzip.min (+ graph.shortread.zipcodes)
  # rather than a plain graph.min.  We check for any file ending in .min.
  local min_file

  mkdir -p "${outdir}"

  min_file=$(ls "${outdir}"/graph*.min 2>/dev/null | head -1 || true)
  if [[ -s "${gbz}" && -s "${dist}" && -n "${min_file}" && -s "${min_file}" ]]; then
    echo "[$(date)] Reusing cached vg index: ${outdir}"
    return 0
  fi

  (
    flock -x 200
    min_file=$(ls "${outdir}"/graph*.min 2>/dev/null | head -1 || true)
    if [[ -s "${gbz}" && -s "${dist}" && -n "${min_file}" && -s "${min_file}" ]]; then
      echo "[$(date)] Reusing cached vg index (built by parallel job): ${outdir}"
      exit 0
    fi
    echo "[$(date)] Building vg giraffe index for DEL ${size}"
    echo "[$(date)]   ref:    ${REF}"
    echo "[$(date)]   vcf:    ${vcf}"
    echo "[$(date)]   prefix: ${outdir}/graph"

    vg autoindex \
      --workflow giraffe \
      -r "${REF}" \
      -v "${vcf}" \
      -p "${outdir}/graph" \
      -t "${SLURM_CPUS_PER_TASK:-4}"

    echo "[$(date)] Index built:"
    ls -lh "${outdir}/"
  ) 200>"${outdir}/.build.lock"
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
case "${SV_TYPE}" in
  "DEL")
    for SIZE in "${!DEL_SIZES[@]}"; do

      VCF="${VCF_DIR}/del_${SIZE}.vcf.gz"
      READS_DIR="${READS}/${COV_LABEL}/var_del_${SIZE}_${RUN_TAG}"
      BASE_DIR="${RESULTS}/${COV_LABEL}_${ERR_LABEL_DIR}/${SIZE}/${N_LABEL}/${FREQ_LABEL_DIR}/${K_LABEL}"
      mkdir -p "${BASE_DIR}"

      PREFIX="var_del_${SIZE}_${RUN_TAG}.vg_giraffe"
      GAM="${BASE_DIR}/${PREFIX}.aln.gam"
      PACK="${BASE_DIR}/${PREFIX}.pack"
      CALLS="${BASE_DIR}/${PREFIX}.calls.vcf"
      AF_OUT="${BASE_DIR}/${PREFIX}.allele_frequencies.k${K}.tsv"
      RUNTIMES="${BASE_DIR}/${PREFIX}.runtimes.k${K}.tsv"

      echo
      echo "==== DEL ${SIZE}  (RUN_TAG=${RUN_TAG}) ===="
      echo "     VCF:      ${VCF}"
      echo "     READS:    ${READS_DIR}"
      echo "     BASE_DIR: ${BASE_DIR}"

      if [[ -s "${AF_OUT}" ]]; then
        echo "[$(date)] Skipping — AF already exists: ${AF_OUT}"
        continue
      fi

      # Preflight
      [[ -s "${REF}" && -s "${REF}.fai" ]] \
        || { echo "ERROR: missing REF or .fai: ${REF}" >&2; exit 1; }
      [[ -s "${VCF}" && -s "${VCF}.tbi" ]] \
        || { echo "SKIP: missing VCF/.tbi for ${SIZE}: ${VCF}" >&2; continue; }
      [[ -s "${READS_DIR}/r1.fq" && -s "${READS_DIR}/r2.fq" ]] \
        || { echo "ERROR: missing reads at ${READS_DIR}" >&2; exit 1; }

      # Per-step runtime TSV (sibling of AF output): method size step elapsed_s cached
      printf "method\tsize\tstep\telapsed_s\tcached\n" > "${RUNTIMES}.tmp"
      mv -f "${RUNTIMES}.tmp" "${RUNTIMES}"

      # Index (cached per pos_label × SIZE — first job pays the build cost)
      INDEX_DIR="${VG_INDEX_ROOT}/${SIZE}"
      _t0=${SECONDS}
      _idx_cached=1
      [[ -s "${INDEX_DIR}/graph.giraffe.gbz" ]] || _idx_cached=0
      ensure_index "${SIZE}" "${VCF}" "${INDEX_DIR}"
      printf "vg_giraffe\t%s\tautoindex\t%d\t%d\n" "${SIZE}" $(( SECONDS - _t0 )) "${_idx_cached}" >> "${RUNTIMES}"

      GBZ="${INDEX_DIR}/graph.giraffe.gbz"
      DIST="${INDEX_DIR}/graph.dist"
      # vg 1.70 may name these graph.shortread.withzip.min / graph.shortread.zipcodes
      MIN=$(ls "${INDEX_DIR}"/graph*.min 2>/dev/null | head -1)
      if [[ -z "${MIN}" || ! -s "${MIN}" ]]; then
        echo "ERROR: no minimizer index (*.min) found in ${INDEX_DIR}" >&2
        exit 1
      fi
      ZIPCODES_ARG=""
      ZIPCODES_FILE=$(ls "${INDEX_DIR}"/graph*.zipcodes 2>/dev/null | head -1 || true)
      if [[ -n "${ZIPCODES_FILE}" && -s "${ZIPCODES_FILE}" ]]; then
        ZIPCODES_ARG="-z ${ZIPCODES_FILE}"
      fi

      # Map (skip if GAM already on disk from a previous run)
      _t0=${SECONDS}
      if [[ -s "${GAM}" ]]; then
        echo "[$(date)] Reusing existing GAM: ${GAM}"
        printf "vg_giraffe\t%s\tgiraffe\t0\t1\n" "${SIZE}" >> "${RUNTIMES}"
      else
        echo "[$(date)] vg giraffe  (threads=${SLURM_CPUS_PER_TASK:-4})"
        vg giraffe \
          -Z "${GBZ}" -d "${DIST}" -m "${MIN}" ${ZIPCODES_ARG} \
          -f "${READS_DIR}/r1.fq" -f "${READS_DIR}/r2.fq" \
          -t "${SLURM_CPUS_PER_TASK:-4}" \
          --progress \
          > "${GAM}.tmp"
        mv -f "${GAM}.tmp" "${GAM}"
        ls -lh "${GAM}"
        printf "vg_giraffe\t%s\tgiraffe\t%d\t0\n" "${SIZE}" $(( SECONDS - _t0 )) >> "${RUNTIMES}"
      fi

      # Pack read support onto graph (skip if already on disk)
      _t0=${SECONDS}
      if [[ -s "${PACK}" ]]; then
        echo "[$(date)] Reusing existing pack: ${PACK}"
        printf "vg_giraffe\t%s\tpack\t0\t1\n" "${SIZE}" >> "${RUNTIMES}"
      else
        echo "[$(date)] vg pack"
        vg pack \
          -x "${GBZ}" \
          -g "${GAM}" \
          -o "${PACK}.tmp" \
          -e \
          -t "${SLURM_CPUS_PER_TASK:-4}"
        mv -f "${PACK}.tmp" "${PACK}"
        ls -lh "${PACK}"
        printf "vg_giraffe\t%s\tpack\t%d\t0\n" "${SIZE}" $(( SECONDS - _t0 )) >> "${RUNTIMES}"
      fi

      # Genotype the graph's snarls from read support.  We intentionally
      # *don't* pass -v here: the graph was built (via autoindex) from a
      # single-variant truth VCF, so vg call emits exactly the deletion
      # snarl we want.  Passing -v requires a haplotype-tracked graph
      # (produced by `vg construct -a`, not autoindex), and silently
      # returns no records when that annotation is missing.
      echo "[$(date)] vg call"
      _t0=${SECONDS}
      vg call "${GBZ}" \
        -k "${PACK}" \
        -t "${SLURM_CPUS_PER_TASK:-4}" \
        > "${CALLS}.tmp"
      mv -f "${CALLS}.tmp" "${CALLS}"
      printf "vg_giraffe\t%s\tcall\t%d\t0\n" "${SIZE}" $(( SECONDS - _t0 )) >> "${RUNTIMES}"
      echo "[$(date)] vg call output:"
      grep -v "^##" "${CALLS}"

      # Extract AF for the deletion site.  The truth VCF has exactly one record,
      # so match on CHROM+POS to be robust to extra sites vg may emit.
      echo "[$(date)] Extracting AF from AD"
      TRUTH_POS=$(bcftools query -f '%POS\n' "${VCF}" | head -1)
      if [[ -z "${TRUTH_POS}" ]]; then
        echo "ERROR: could not read POS from truth VCF ${VCF}" >&2
        exit 1
      fi

      # The graph was built from a single-variant VCF, so vg call normally
      # emits exactly one record.  Prefer the record whose POS matches the
      # truth VCF; if there is no such match (autoindex may nudge the
      # anchor), fall back to the first data record.
      AF_LINE=$(bcftools query -f '%CHROM\t%POS\t[%AD]\n' "${CALLS}" \
                | awk -v chrom="${CHROM}" -v pos="${TRUTH_POS}" '
                    function emit(ad_str,    n, ad, ref, alt, i, tot) {
                      n = split(ad_str, ad, ",")
                      ref = ad[1] + 0
                      alt = 0
                      for (i = 2; i <= n; i++) alt += ad[i] + 0
                      tot = ref + alt
                      if (tot == 0) { printf "NA|NA"; return }
                      printf "%.7f|%.7f", 1 - alt/tot, alt/tot
                    }
                    $1==chrom && $2==pos { emit($3); matched=1; exit }
                    # Remember first data row in case no POS match is found
                    NR==1 { first_ad=$3 }
                    END {
                      if (matched) ;
                      else if (first_ad != "") emit(first_ad)
                      else printf "NA|NA"
                    }')

      if [[ -z "${AF_LINE}" ]]; then
        echo "ERROR: empty AF extraction" >&2
        exit 1
      fi

      printf "%s\n" "${AF_LINE}" > "${AF_OUT}.tmp"
      mv -f "${AF_OUT}.tmp" "${AF_OUT}"

      echo
      echo "[$(date)] Done for DEL ${SIZE}.  AF (af_ref|af_alt):"
      cat "${AF_OUT}"
    done
    ;;
  *)
    echo "05b_vg_giraffe.sh: unsupported SV_TYPE=${SV_TYPE} (expected DEL)" >&2
    exit 1
    ;;
esac
