#!/bin/bash

# =============================================================================
# config_sv_var_deletions.sh
# Configuration for the variation-aware pipeline:
#   N ecotypes from GrENET VCF → per-sample consensus FASTAs →
#   N_SV of them get a deletion → SHORtS pools all N at equal fractions
# =============================================================================

# Common paths
WORK=/home/tbellagio/scratch/visor_freqk
REF=${WORK}/data/reference/Chr1.fa

# GrENET VCF (SNPs only, samples named as ecotype IDs)
# Original VCF (chrom = "1") — run 00_prep_vcf.sh once to produce RENAMED_VCF
VCF_FILE=${WORK}/data/reference/greneNet_final_v1.1.recode.vcf.gz

# Chromosome name mapping file: VCF "1" → FASTA "Chr1"
CHROM_MAP=${WORK}/data/reference/chrom_map.txt

# Renamed VCF (chrom = "Chr1") — produced by 00_prep_vcf.sh
# All downstream scripts use this, not VCF_FILE directly
RENAMED_VCF=${WORK}/data/reference/greneNet_final_v1.1.recode.Chr1.vcf.gz

# SV design: deletions on Chr1 around 10 Mb
SV_TYPE="DEL"
CHROM="Chr1"
SV_START_0=10000000      # BED 0-based start
ANCHOR_POS=9999999       # one base before deletion (VCF POS, 1-based)
POS_LABEL="pos10mb"

# Test: use first N_SAMPLES ecotypes from the VCF header
N_SAMPLES=10

# Fraction of the N_SAMPLES haplotypes that carry the SV
# N_SV = round(N_SAMPLES * SV_FREQ); must be >=1 and < N_SAMPLES
SV_FREQ=0.50

# Only one SV size for the variation test
declare -A DEL_SIZES=(
  ["1kb"]=1000
)

# Pool-seq design
COVERAGE=50
ERROR_RATE=0.001
K=31
# FREQ is an alias for SV_FREQ so that 04_make_vcf.sh / 05_freqk.sh
# (which use ${FREQ}) work unchanged with this config
FREQ=${SV_FREQ}

# ---------------------------------------------------------------------------
# Directory roots
#
#   BEDS      — shared with old pipeline (BED content is position-only, identical)
#   HAPS_VAR  — new: per-ecotype consensus FASTAs + SV-injected FASTAs
#               (old pipeline uses  data/haplotypes/ — no overlap)
#   READS_VAR — new reads dir
#               (old pipeline uses  data/reads/      — no overlap)
#   READS     — alias for READS_VAR so that 05_freqk.sh (which uses ${READS})
#               finds the right read files
#   VCF_DIR   — /var subdirectory: keeps SV VCFs separate from old pipeline
#               (future runs may produce different VCFs, e.g. INS, other positions)
#   RESULTS   — /var subdirectory: freqk sees different reads → different outputs
# ---------------------------------------------------------------------------
HAPS_VAR=${WORK}/data/haplotypes_var/del/${POS_LABEL}
BEDS=${WORK}/data/beds/del/${POS_LABEL}
READS_VAR=${WORK}/data/reads_var/del/${POS_LABEL}
READS=${READS_VAR}                                           # alias for 05_freqk.sh
VCF_DIR=${WORK}/data/vcf/del/${POS_LABEL}/var
RESULTS=${WORK}/results/del/${POS_LABEL}/var

FREQK=/home/tbellagio/scratch/bin/freqk
