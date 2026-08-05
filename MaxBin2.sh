#!/bin/bash
# Usage:
#   bash MaxBin2.sh <SampleName>

SAMPLE="$1"

CONTIGS="${SAMPLE}_assembly/metaspades/renamed_${SAMPLE}_1k_contigs.fasta"
R1="${SAMPLE}_host_unmapped/${SAMPLE}_paired_unmapped_R1.fq.gz"
R2="${SAMPLE}_host_unmapped/${SAMPLE}_paired_unmapped_R2.fq.gz"
U="${SAMPLE}_host_unmapped/${SAMPLE}_unpaired_unmapped.fq.gz"

OUTDIR="${SAMPLE}_maxbin2"
mkdir -p "${OUTDIR}"

cat "${R1}" "${R2}" "${U}" > "${OUTDIR}/merged.fq.gz"

run_MaxBin.pl \
  -contig "${CONTIGS}" \
  -reads "${OUTDIR}/merged.fq.gz" \
  -out   "${OUTDIR}/output"
