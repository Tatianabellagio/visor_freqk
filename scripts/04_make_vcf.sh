#!/bin/bash
#SBATCH --job-name=make_vcf_all
#SBATCH --output=logs/04_make_vcf_%j.out
#SBATCH --error=logs/04_make_vcf_%j.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
source "$(mamba info --base)/etc/profile.d/conda.sh" && conda activate pang

# =============================================================================
# 04_make_vcf.sh
# Purpose: Create VCFs for all SV sizes defined in the config
#          - DEL: sequence-resolved VCFs (multi-size DELs around Chr1:10Mb),
#                  using the same anchor logic as the validated 1kb case.
#          - INS: simple symbolic insertion VCFs (ALT=<INS>, SVLEN from INS_SIZES).
# =============================================================================

CONFIG_FILE=${1:-"$(dirname "$0")/config_sv_deletions.sh"}
source "${CONFIG_FILE}"

# VCF output is per (SV_TYPE, SIZE) only; content does not depend on COVERAGE, FREQ, or ERROR_RATE.
# Different experiments (cov/err/freq) share the same VCF; only reads and results paths vary.
mkdir -p "${VCF_DIR}" logs

case "${SV_TYPE}" in
  "DEL")
    for SIZE in "${!DEL_SIZES[@]}"; do
      DEL_LEN=${DEL_SIZES[$SIZE]}

      # Calculate coordinates
      # BED: CHROM SV_START_0 SV_START_0+DEL_LEN (0-based half-open)
      # VCF: POS=ANCHOR_POS (anchor), DEL_END=SV_START_0+DEL_LEN (1-based)
      POS=${ANCHOR_POS}
      DEL_END=$((SV_START_0 + DEL_LEN))
      SVLEN=$((DEL_LEN + 1))  # anchor base + deleted bases

      VCF_RAW="${VCF_DIR}/del_${SIZE}.vcf"
      VCF_GZ="${VCF_DIR}/del_${SIZE}.vcf.gz"

      echo "[$(date)] Building DEL VCF for ${SIZE}: POS=${POS}, DEL_END=${DEL_END}, SVLEN=${SVLEN}"

      # Extract REF sequence (anchor + deleted region)
      REF_SEQ=$(samtools faidx "${REF}" "${CHROM}:${POS}-${DEL_END}" | grep -v "^>" | tr -d '\n')
      ALT_SEQ="${REF_SEQ:0:1}"  # anchor base only

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

      # For insertions, use a sequence-resolved representation:
      # POS:  1-based coordinate of insertion site (same anchor as in BED)
      # REF:  reference base at POS
      # ALT:  REF base + inserted sequence (must match BED's inserted sequence)
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
    echo "04_make_vcf.sh: unsupported SV_TYPE=${SV_TYPE} (expected DEL or INS)" >&2
    exit 1
    ;;
esac

echo "[$(date)] All VCFs created successfully"