#!/bin/bash

# Usage:
#   bash trim.sh SampleName
#

file1="${1}_1.fq.gz"
file2="${1}_2.fq.gz"

outdir="trimmed"
mkdir -p "${outdir}"

java -jar /path/to/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 1 -phred33 \
  "${file1}" "${file2}" \
  "${outdir}/paired_1.fq" "${outdir}/unpaired_1.fq" \
  "${outdir}/paired_2.fq" "${outdir}/unpaired_2.fq" \
  ILLUMINACLIP:/path/to/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 \
  SLIDINGWINDOW:10:15 MINLEN:50

cat "${outdir}/unpaired_1.fq" "${outdir}/unpaired_2.fq" > "${outdir}/unpaired.fq"
pigz "${outdir}"/*.fq





