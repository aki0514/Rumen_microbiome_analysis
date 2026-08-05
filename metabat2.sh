#!/bin/bash
# Usage:
#   bash metabat2.sh <SampleName>
#

SAMPLE="$1"

CONTIGS="${SAMPLE}_assembly/metaspades/renamed_${SAMPLE}_1k_contigs.fasta"
BAM="${SAMPLE}_backmap/${SAMPLE}_BWA_back.bam"
OUTDIR="${SAMPLE}_metabat2"
mkdir -p "${OUTDIR}"

jgi_summarize_bam_contig_depths \
  --outputDepth "${OUTDIR}/${SAMPLE}_depth.txt" \
  "${BAM}"

metabat2 \
  -i "${CONTIGS}" \
  -a "${OUTDIR}/${SAMPLE}_depth.txt" \
  -o "${OUTDIR}/bin"


