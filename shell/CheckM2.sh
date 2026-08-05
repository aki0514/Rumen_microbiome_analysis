#!/bin/bash
# Usage:
#   bash CheckM2.sh <BinDir> <OutputPrefix>
#
# <BinDir> should contain all DAS Tool refined bins (*.fa) from all samples.


BINS="$1"
PREFIX="$2"

OUTDIR="${PREFIX}_CheckM2"
mkdir -p "${OUTDIR}"

checkm2 predict \
  --input "${BINS}" \
  --output-directory "${OUTDIR}" \
  -x fa

