#!/bin/bash

# Usage:
#   bash coverm_protein.sh <ReferenceFNA> <ForwardReads> <ReverseReads> <OutputDir> <SampleName>
#
# <ReferenceFNA> : non-redundant nucleotide FASTA corresponding to the protein catalog
# <ForwardReads> : forward paired-end reads (FASTQ.gz)
# <ReverseReads> : reverse paired-end reads (FASTQ.gz)
# <OutputDir>    : output directory
# <SampleName>   : sample name

REFERENCE="$1"
READ1="$2"
READ2="$3"
OUTDIR="$4"
SAMPLE="$5"

mkdir -p "${OUTDIR}"

coverm contig \
    -1 "${READ1}" \
    -2 "${READ2}" \
    --reference "${REFERENCE}" \
    -m rpkm \
    -p bwa-mem \
    --threads 16 \
    --min-covered-fraction 80 \
    --min-read-percent-identity 95 \
    --min-read-aligned-percent 90 \
    > "${OUTDIR}/${SAMPLE}_protein.tsv"
