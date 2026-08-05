#!/bin/bash
# Usage:
#   bash dRep99.sh <BinDir> <OutputPrefix> <CheckM2_report>
#
# <BinDir>         : directory containing all bins (*.fa)
# <OutputPrefix>   : prefix for output directory
# <CheckM2_report> : CheckM2 quality report (CSV)

BINS="$1"
PREFIX="$2"
CHECKM2="$3"

OUTDIR="${PREFIX}_dRep"
mkdir -p "${OUTDIR}"

dRep dereplicate "${OUTDIR}" \
  -g "${BINS}"/*.fa \
  -sa 0.99 \
  -comp 50 \
  -con 10 \
  --genomeInfo "${CHECKM2}"
