#!/bin/bash

# Common paths
WORK=/home/tbellagio/scratch/pang/visor_freqk
REF=${WORK}/data/reference/Chr1.fa

# SV design: deletions on Chr1 around 10 Mb
SV_TYPE="DEL"
CHROM="Chr1"
SV_START_0=10000000      # BED 0-based start
ANCHOR_POS=9999999       # one base before deletion (VCF POS, 1-based)
POS_LABEL="pos10mb"      # human-readable label for this position

# Design-dependent roots (per SV_TYPE and POS_LABEL)
BEDS=${WORK}/data/beds/del/${POS_LABEL}
HAPS=${WORK}/data/haplotypes/del/${POS_LABEL}
READS=${WORK}/data/reads/del/${POS_LABEL}
VCF_DIR=${WORK}/data/vcf/del/${POS_LABEL}
RESULTS=${WORK}/results/del/${POS_LABEL}

declare -A DEL_SIZES=(
  ["100bp"]=100
  ["500bp"]=500
  ["1kb"]=1000
  ["5kb"]=5000
  ["10kb"]=10000
)

# Pool-seq design
FREQ=0.50
COVERAGE=50
ERROR_RATE=0.000
# Shared WT clone (reference only, does not depend on POS_LABEL)
WT_CLONE=${WORK}/data/haplotypes/del/_clone_WT

# freqk
K=31
FREQK=/home/tbellagio/scratch/pang/test_freqk/freqk/target/release/freqk

