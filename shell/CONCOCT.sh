#!/bin/bash
# Usage:
#   bash CONCOCT.sh <SampleName>

SAMPLE="$1"

CONTIGS="${SAMPLE}_assembly/metaspades/renamed_${SAMPLE}_1k_contigs.fasta"
BAM="${SAMPLE}_backmap/${SAMPLE}_BWA_back.bam"

OUTDIR="${SAMPLE}_concoct"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

cut_up_fasta.py "${CONTIGS}" -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa

concoct_coverage_table.py contigs_10K.bed "${BAM}" > coverage_table.tsv

concoct --composition_file contigs_10K.fa \
        --coverage_file    coverage_table.tsv \
        -b concoct_output/

merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv
mkdir -p concoct_output/fasta_bins
extract_fasta_bins.py "${CONTIGS}" concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins
