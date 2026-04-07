#!/bin/bash
#SBATCH --job-name=make_vcf_var
#SBATCH --output=logs/04_make_vcf_var_%j.out
#SBATCH --error=logs/04_make_vcf_var_%j.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
# =============================================================================
# 04_make_vcf_var.sh
# Purpose: Create VCFs for all SV sizes defined in the config.
#          Writes to VCF_DIR = data/vcf/del/<pos>/var/ (var-pipeline subdir).
#          Logic is identical to 04_make_vcf.sh; kept as a standalone script
#          so both pipelines can evolve independently.
# =============================================================================
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CONFIG_FILE=${1:-"$(dirname "$0")/config_sv_var_deletions.sh"}
source "${CONFIG_FILE}"

mkdir -p "${VCF_DIR}" logs

# ---------------------------------------------------------------------------
# Fast early-exit: if all VCFs already exist, nothing to do.
# Runs before conda activation to avoid unnecessary overhead.
# ---------------------------------------------------------------------------
_N_VCF_DONE=0; _N_VCF_TOTAL=0
case "${SV_TYPE}" in
  "DEL") for _SZ in "${!DEL_SIZES[@]}"; do
           _N_VCF_TOTAL=$((_N_VCF_TOTAL+1))
           [[ -s "${VCF_DIR}/del_${_SZ}.vcf.gz" && -s "${VCF_DIR}/del_${_SZ}.vcf.gz.tbi" ]] \
             && _N_VCF_DONE=$((_N_VCF_DONE+1))
         done ;;
  "INS") for _SZ in "${!INS_SIZES[@]}"; do
           _N_VCF_TOTAL=$((_N_VCF_TOTAL+1))
           [[ -s "${VCF_DIR}/ins_${_SZ}.vcf.gz" && -s "${VCF_DIR}/ins_${_SZ}.vcf.gz.tbi" ]] \
             && _N_VCF_DONE=$((_N_VCF_DONE+1))
         done ;;
esac
if [[ "${_N_VCF_TOTAL}" -gt 0 && "${_N_VCF_DONE}" -ge "${_N_VCF_TOTAL}" ]]; then
  echo "[$(date)] Skipping 04_make_vcf_var — all ${_N_VCF_TOTAL} VCFs already exist in ${VCF_DIR}"
  exit 0
fi
echo "[$(date)] Found ${_N_VCF_DONE}/${_N_VCF_TOTAL} VCFs — building missing ones."

source "$(mamba info --base)/etc/profile.d/conda.sh" && conda activate pang

case "${SV_TYPE}" in
  "DEL")
    for SIZE in "${!DEL_SIZES[@]}"; do
      DEL_LEN=${DEL_SIZES[$SIZE]}

      POS=${ANCHOR_POS}
      DEL_END=$((SV_START_0 + DEL_LEN))
      SVLEN=$((DEL_LEN + 1))

      VCF_RAW="${VCF_DIR}/del_${SIZE}.vcf"
      VCF_GZ="${VCF_DIR}/del_${SIZE}.vcf.gz"

      # Fast path (no lock)
      if [[ -s "${VCF_GZ}" && -s "${VCF_GZ}.tbi" ]]; then
        echo "[$(date)] Reusing existing VCF for DEL ${SIZE}: ${VCF_GZ}"
        continue
      fi

      # Slow path: acquire per-size lock so concurrent sv_freq jobs don't race
      (
        flock -x 200
        if [[ -s "${VCF_GZ}" && -s "${VCF_GZ}.tbi" ]]; then
          echo "[$(date)] Reusing existing VCF for DEL ${SIZE} (built by parallel job): ${VCF_GZ}"
          exit 0
        fi

      echo "[$(date)] Building DEL VCF for ${SIZE}: POS=${POS}, DEL_END=${DEL_END}, SVLEN=${SVLEN}"

      REF_SEQ=$(samtools faidx "${REF}" "${CHROM}:${POS}-${DEL_END}" | grep -v "^>" | tr -d '\n')
      ALT_SEQ="${REF_SEQ:0:1}"

      echo "[$(date)]   REF length: ${#REF_SEQ}, ALT: ${ALT_SEQ}"

      cat > "${VCF_RAW}" << VCFEOF
##fileformat=VCFv4.2
##source=simulated_visor
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##contig=<ID=${CHROM}>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
VCFEOF

      echo -e "${CHROM}\t${POS}\t.\t${REF_SEQ}\t${ALT_SEQ}\t1000\tPASS\tSVTYPE=DEL;SVLEN=-${SVLEN};END=${DEL_END}" >> "${VCF_RAW}"

      bgzip -f "${VCF_RAW}"
      tabix "${VCF_GZ}"

      echo "[$(date)] Done: ${VCF_GZ}"
      ) 200>"${VCF_DIR}/del_${SIZE}.lock"
    done
    ;;
  "INS")
    for SIZE in "${!INS_SIZES[@]}"; do
      INS_LEN=${INS_SIZES[$SIZE]}

      POS=${SV_START_0}
      END=${POS}
      SVLEN=${INS_LEN}

      VCF_RAW="${VCF_DIR}/ins_${SIZE}.vcf"
      VCF_GZ="${VCF_DIR}/ins_${SIZE}.vcf.gz"

      # Fast path (no lock)
      if [[ -s "${VCF_GZ}" && -s "${VCF_GZ}.tbi" ]]; then
        echo "[$(date)] Reusing existing VCF for INS ${SIZE}: ${VCF_GZ}"
        continue
      fi

      # Slow path: acquire per-size lock
      (
        flock -x 200
        if [[ -s "${VCF_GZ}" && -s "${VCF_GZ}.tbi" ]]; then
          echo "[$(date)] Reusing existing VCF for INS ${SIZE} (built by parallel job): ${VCF_GZ}"
          exit 0
        fi

      echo "[$(date)] Building INS VCF for ${SIZE}: POS=${POS}, SVLEN=${SVLEN}"

      REF_BASE=$(samtools faidx "${REF}" "${CHROM}:${POS}-${POS}" | grep -v "^>" | tr -d '\n')
      BED="${BEDS}/hack_ins_${SIZE}.bed"
      INS_SEQ=$(awk 'NR==1 {print $5}' "${BED}")
      ALT_SEQ="${REF_BASE}${INS_SEQ}"

      cat > "${VCF_RAW}" << VCFEOF
##fileformat=VCFv4.2
##source=simulated_visor
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##contig=<ID=${CHROM}>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
VCFEOF

      echo -e "${CHROM}\t${POS}\t.\t${REF_BASE}\t${ALT_SEQ}\t1000\tPASS\tSVTYPE=INS;SVLEN=${SVLEN};END=${END}" >> "${VCF_RAW}"

      bgzip -f "${VCF_RAW}"
      tabix "${VCF_GZ}"

      echo "[$(date)] Done: ${VCF_GZ}"
      ) 200>"${VCF_DIR}/ins_${SIZE}.lock"
    done
    ;;
  *)
    echo "04_make_vcf_var.sh: unsupported SV_TYPE=${SV_TYPE} (expected DEL or INS)" >&2
    exit 1
    ;;
esac

echo "[$(date)] All VCFs created successfully"
