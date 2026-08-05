#!/bin/bash
# Usage:
#   bash Assemble_MAG.sh <SampleName>

SAMPLE="$1"

INDIR="${SAMPLE}_host_unmapped"
R1="${INDIR}/${SAMPLE}_paired_unmapped_R1.fq.gz"
R2="${INDIR}/${SAMPLE}_paired_unmapped_R2.fq.gz"
U="${INDIR}/${SAMPLE}_unpaired_unmapped.fq.gz"

OUTDIR="${SAMPLE}_assembly"
mkdir -p "${OUTDIR}"

spades.py \
  -1 "${R1}" \
  -2 "${R2}" \
  -s "${U}" \
  -o "${OUTDIR}/metaspades" \
  -m 200 --meta -k auto

seqkit seq -m 1000 \
  "${OUTDIR}/metaspades/contigs.fasta" \
  > "${OUTDIR}/metaspades/${SAMPLE}_1k_contigs.fasta"

awk -v sample="${SAMPLE}" '/^>/{sub(/^>NODE/, ">"sample, $0); print; next}{print}' \
  "${OUTDIR}/metaspades/${SAMPLE}_1k_contigs.fasta" \
  > "${OUTDIR}/metaspades/renamed_${SAMPLE}_1k_contigs.fasta"
