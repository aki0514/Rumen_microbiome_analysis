#!/bin/bash
# Usage:
#   bash DAS_Tool.sh <SampleName>

SAMPLE="$1"

# paths
CONTIGS="${SAMPLE}_assembly/metaspades/renamed_${SAMPLE}_1k_contigs.fasta"
CONCOCT_DIR="${SAMPLE}_concoct/concoct_output/fasta_bins"
MAXBIN_DIR="${SAMPLE}_maxbin2"
METABAT_DIR="${SAMPLE}_metabat2"

OUTDIR="${SAMPLE}_DASTool"
mkdir -p "${OUTDIR}"

Fasta_to_Contig2Bin.sh -i "${CONCOCT_DIR}" -e fa     > "${OUTDIR}/concoct.tsv"
Fasta_to_Contig2Bin.sh -i "${MAXBIN_DIR}"  -e fasta  > "${OUTDIR}/maxbin2.tsv"
Fasta_to_Contig2Bin.sh -i "${METABAT_DIR}" -e fa     > "${OUTDIR}/metabat2.tsv"

DAS_Tool \
  -i "${OUTDIR}/concoct.tsv,${OUTDIR}/maxbin2.tsv,${OUTDIR}/metabat2.tsv" \
  -l concoct,maxbin2,metabat2 \
  -c "${CONTIGS}" \
  -o "${OUTDIR}/${SAMPLE}" \
  --write_bins
