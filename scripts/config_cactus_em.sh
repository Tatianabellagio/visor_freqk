#!/bin/bash
# =============================================================================
# config_cactus_em.sh — VISOR sim with cactus founder FASTAs (no GrENE-VCF
# consensus step) + cactus pangenome k-mer EM pipeline (no freqk).
#
# Differences from config_sv_var_deletions.sh:
#   - WT haplotypes come from cactus founder FASTAs directly (not bcftools
#     consensus). Skip 00_prep_vcf and 00_apply_vcf.
#   - Sample IDs are cactus accession IDs (e.g. "100042"), set from the
#     cactus FASTA filenames.
#   - Reads downstream go to our cactus-em pipeline (05c_cactus_em.sh)
#     instead of (or alongside) freqk (05_freqk_var.sh).
# =============================================================================

# Paths
WORK=/home/tbellagio/scratch/visor_freqk          # leverage existing visor_freqk infra
CACTUS_DIR=/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only
REF=${CACTUS_DIR}/TAIR10.chr.fa

# Pool composition
N_SAMPLES=10
SV_FREQ=0.5

# SV design
SV_TYPE="DEL"
CHROM="Chr1"
SV_START_0=10000000
ANCHOR_POS=9999999
POS_LABEL="cactus_pos10mb"

declare -A DEL_SIZES=(
  ["1kb"]=1000
)

# Pool-seq design
COVERAGE=30
ERROR_RATE=0.001
K=31
READ_LEN=150
FREQ=${SV_FREQ}

# Directory roots (separate /cactus_em subdir to keep clean from freqk runs)
HAPS_WT=${WORK}/data/haplotypes_var/cactus_em/${POS_LABEL}/wt_n${N_SAMPLES}
HAPS_SV=${WORK}/data/haplotypes_var/cactus_em/${POS_LABEL}
BEDS=${WORK}/data/beds/cactus_em/${POS_LABEL}
READS_VAR=${WORK}/data/reads_var/cactus_em/${POS_LABEL}
READS=${READS_VAR}
VCF_DIR=${WORK}/data/vcf/cactus_em/${POS_LABEL}/var
RESULTS=${WORK}/results/cactus_em/${POS_LABEL}/var

# Our cactus+EM pipeline
POOLFREQ=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq
CN_KMER_PREFIX=${POOLFREQ}/data/cn_full
CN_VAR=${POOLFREQ}/data/cn_var_82.cn_var.npz
CN_VAR_META=${POOLFREQ}/data/cn_var_82.meta.npz
PYTHON_HAPFM=/home/tbellagio/miniforge3/envs/hapfm/bin/python

# Region BED for SHORtS — Chr1 only (matches visor_freqk default scope)
REGION_BED_PATH=${BEDS}/region_chr1.bed
