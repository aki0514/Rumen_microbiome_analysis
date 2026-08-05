q#!/bin/bash
# Usage:
#   bash prodigal.sh <GenomeDir> <OutputPrefix>
#
# <GenomeDir>   : directory containing dereplicated genomes (*.fa)
# <OutputPrefix>: prefix for output directory

GENOMES="$1"
PREFIX="$2"

OUTDIR="${PREFIX}_prodigal"
mkdir -p "${OUTDIR}"

for FILE in "${GENOMES}"/*.fa; do
    BASENAME=$(basename "${FILE}" .fa)

    prodigal \
        -f gff \
        -i "${FILE}" \
        -o "${OUTDIR}/${BASENAME}.gff" \
        -a "${OUTDIR}/${BASENAME}_protein.faa" \
        -d "${OUTDIR}/${BASENAME}_nucleotide.fna" \
        -p meta

    # Filter proteins with complete ORFs (remove partial=xx except 00)
    awk '/^>/{p=($0 ~ /partial=00/)} p' \
        "${OUTDIR}/${BASENAME}_protein.faa" \
        > "${OUTDIR}/${BASENAME}_protein_filtered.faa"
done
