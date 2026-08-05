#!/bin/bash
# Usage:
#   bash BWA_reference.sh <Reference.fa> <SampleName> <R1.fq.gz> <R2.fq.gz> <U.fq.gz>
#
# Example:
#   bash BWA_reference.sh cow_ref.fa SampleA \
#       trimmed/SampleA_paired_output_1.fq.gz \
#       trimmed/SampleA_paired_output_2.fq.gz \
#       trimmed/SampleA_unpaired_output.fq.gz


REF="$1"
SAMPLE="$2"
R1="$3"
R2="$4"
U="$5"

OUT="${SAMPLE}_host_unmapped"
mkdir -p "${OUT}"

# --- BWA mapping ---
bwa mem "${REF}" "${R1}" "${R2}" > "${OUT}/${SAMPLE}_paired.sam"
bwa mem "${REF}" "${U}" > "${OUT}/${SAMPLE}_unpaired.sam"

samtools sort -@ 8 -O bam -o "${OUT}/${SAMPLE}_paired.bam"   "${OUT}/${SAMPLE}_paired.sam"
samtools sort -@ 8 -O bam -o "${OUT}/${SAMPLE}_unpaired.bam" "${OUT}/${SAMPLE}_unpaired.sam"

samtools index "${OUT}/${SAMPLE}_paired.bam"
samtools index "${OUT}/${SAMPLE}_unpaired.bam"

samtools view -b -f 4 -F 0x800 "${OUT}/${SAMPLE}_paired.bam"   > "${OUT}/${SAMPLE}_paired_unmapped.bam"
samtools view -b -f 4 -F 0x800 "${OUT}/${SAMPLE}_unpaired.bam" > "${OUT}/${SAMPLE}_unpaired_unmapped.bam"

samtools fastq -s "${OUT}/single.fq" \
  -1 "${OUT}/${SAMPLE}_paired_unmapped_R1.fq" \
  -2 "${OUT}/${SAMPLE}_paired_unmapped_R2.fq" \
  "${OUT}/${SAMPLE}_paired_unmapped.bam"

samtools fastq -0 "${OUT}/${SAMPLE}_unpaired_unmapped.fq" \
  "${OUT}/${SAMPLE}_unpaired_unmapped.bam"

cat "${OUT}/single.fq" "${OUT}/${SAMPLE}_unpaired_unmapped.fq" > "${OUT}/${SAMPLE}_unpaired_unmapped2.fq"
rm  "${OUT}/single.fq"
rm  "${OUT}/${SAMPLE}_unpaired_unmapped.fq"
mv  "${OUT}/${SAMPLE}_unpaired_unmapped2.fq" "${OUT}/${SAMPLE}_unpaired_unmapped.fq"

pigz "${OUT}"/*.fq

