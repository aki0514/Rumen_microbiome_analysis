#!/bin/bash

# Usage:
#   bash dRep95.sh <BinDir> <OutputPrefix>
#
# <BinDir>       : directory containing genome bins (*.fa)
# <OutputPrefix> : prefix for output directory

BINS="$1"
PREFIX="$2"

OUTDIR="${PREFIX}_dRep95"
mkdir -p "${OUTDIR}"

dRep dereplicate "${OUTDIR}" \
    -g "${BINS}"/*.fa \
    -sa 0.95 \
    --ignoreGenomeQuality \
    -p 8
