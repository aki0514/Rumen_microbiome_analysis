#!/bin/bash

# Usage:
#   bash kofamscan.sh <ProteinFAA> <OutputDir> <Prefix>
#
# <ProteinFAA> : non-redundant protein FASTA file
# <OutputDir>  : output directory
# <Prefix>     : output prefix

PROTEIN="$1"
OUTDIR="$2"
PREFIX="$3"

mkdir -p "${OUTDIR}"

# Run KofamScan
exec_annotation \
    -o "${OUTDIR}/${PREFIX}_kofamscan.txt" \
    "${PROTEIN}" \
    -f mapper
