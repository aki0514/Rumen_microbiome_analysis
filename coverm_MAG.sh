#!/bin/bash

# Usage:
#   bash coverm_MAG.sh <GenomeDir> <ForwardReads> <ReverseReads> <OutputDir> <SampleName>
#
# <GenomeDir>    : directory containing dereplicated MAGs (*.fa)
# <ForwardReads> : forward paired-end reads (FASTQ.gz)
# <ReverseReads> : reverse paired-end reads (FASTQ.gz)
# <OutputDir>    : output directory
# <SampleName>   : sample name

GENOMES="$1"
READ1="$2"
READ2="$3"
OUTDIR="$4"
SAMPLE="$5"

mkdir -p "${OUTDIR}"

coverm genome \
    -1 "${READ1}" \
    -2 "${READ2}" \
    --genome-fasta-directory "${GENOMES}" \
    -x fa \
    -m rpkm \
    -p bwa-mem \
    --threads 16 \
    --min-covered-fraction 10 \
    --min-read-percent-identity 95 \
    --min-read-aligned-percent 75 \
    > "${OUTDIR}/${SAMPLE}_coverm.tsv"
