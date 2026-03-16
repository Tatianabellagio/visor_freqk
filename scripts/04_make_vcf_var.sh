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
source "$(mamba info --base)/etc/profile.d/conda.sh" && conda activate pang

CONFIG_FILE=${1:-"$(dirname "$0")/config_sv_var_deletions.sh"}
source "${CONFIG_FILE}"

mkdir -p "${VCF_DIR}" logs

case "${SV_TYPE}" in
  "DEL")
    for SIZE in "${!DEL_SIZES[@]}"; do
      DEL_LEN=${DEL_SIZES[$SIZE]}

      POS=${ANCHOR_POS}
      DEL_END=$((SV_START_0 + DEL_LEN))
      SVLEN=$((DEL_LEN + 1))

      VCF_RAW="${VCF_DIR}/del_${SIZE}.vcf"
      VCF_GZ="${VCF_DIR}/del_${SIZE}.vcf.gz"

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
    done
    ;;
  *)
    echo "04_make_vcf_var.sh: unsupported SV_TYPE=${SV_TYPE} (expected DEL or INS)" >&2
    exit 1
    ;;
esac

echo "[$(date)] All VCFs created successfully"
