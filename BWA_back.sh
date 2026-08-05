#!/bin/bash
# Usage:
#   bash BWA_back.sh <SampleName>
#

SAMPLE="$1"
CONTIGS="${SAMPLE}_assembly/metaspades/renamed_${SAMPLE}_1k_contigs.fasta"
R1="${SAMPLE}_host_unmapped/${SAMPLE}_paired_unmapped_R1.fq.gz"
R2="${SAMPLE}_host_unmapped/${SAMPLE}_paired_unmapped_R2.fq.gz"
U="${SAMPLE}_host_unmapped/${SAMPLE}_unpaired_unmapped.fq.gz"

OUTDIR="${SAMPLE}_backmap"
mkdir -p "${OUTDIR}"

bwa mem "${CONTIGS}" "${R1}" "${R2}" > "${OUTDIR}/${SAMPLE}_BWA_back_paired.sam"
bwa mem "${CONTIGS}" "${U}"              > "${OUTDIR}/${SAMPLE}_BWA_back_unpaired.sam"

samtools sort -@ 8 -O bam -o "${OUTDIR}/${SAMPLE}_BWA_back_paired.bam"   "${OUTDIR}/${SAMPLE}_BWA_back_paired.sam"
samtools sort -@ 8 -O bam -o "${OUTDIR}/${SAMPLE}_BWA_back_unpaired.bam" "${OUTDIR}/${SAMPLE}_BWA_back_unpaired.sam"

samtools merge -@ 8 "${OUTDIR}/${SAMPLE}_BWA_back.merged.bam" \
  "${OUTDIR}/${SAMPLE}_BWA_back_paired.bam" \
  "${OUTDIR}/${SAMPLE}_BWA_back_unpaired.bam"

samtools view -b -F 0x800 "${OUTDIR}/${SAMPLE}_BWA_back.merged.bam" > "${OUTDIR}/${SAMPLE}_BWA_back.tmp.bam"
samtools index "${OUTDIR}/${SAMPLE}_BWA_back.tmp.bam"

mv "${OUTDIR}/${SAMPLE}_BWA_back.tmp.bam" "${OUTDIR}/${SAMPLE}_BWA_back.bam"


