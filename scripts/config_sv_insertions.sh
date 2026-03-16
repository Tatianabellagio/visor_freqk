#!/bin/bash

# Example config for insertion experiments

WORK=/home/tbellagio/scratch/pang/visor_freqk
REF=${WORK}/data/reference/Chr1.fa

SV_TYPE="INS"
CHROM="Chr1"

# For insertions you'd typically specify the insertion site (0-based) and lengths
SV_START_0=10000000
POS_LABEL="pos10mb"      # human-readable label for this position

# Design-dependent roots (per SV_TYPE and POS_LABEL)
BEDS=${WORK}/data/beds/ins/${POS_LABEL}
HAPS=${WORK}/data/haplotypes/ins/${POS_LABEL}
READS=${WORK}/data/reads/ins/${POS_LABEL}
VCF_DIR=${WORK}/data/vcf/ins/${POS_LABEL}
RESULTS=${WORK}/results/ins/${POS_LABEL}

declare -A INS_SIZES=(
  ["100bp"]=100
  ["500bp"]=500
  ["1kb"]=1000
  ["5kb"]=5000
  ["10kb"]=10000
)

FREQ=0.50
COVERAGE=50
ERROR_RATE=0.000
# Shared WT clone (reference only, does not depend on POS_LABEL)
WT_CLONE=${WORK}/data/haplotypes/ins/_clone_WT

K=31
FREQK=/home/tbellagio/scratch/pang/test_freqk/freqk/target/release/freqk

