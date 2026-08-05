#!/bin/bash
# Usage:
#   bash HMM_cazy.sh <ProteinDir> <OutputPrefix> <dbCAN_HMM> <ParserScript>
#
# <ProteinDir>   : directory containing *_protein_filtered.faa files
# <OutputPrefix> : prefix for output directory
# <dbCAN_HMM>    : path to dbCAN HMM database (e.g., dbCAN-fam-HMMs.txt.v14)
# <ParserScript> : path to hmmscan-parser.sh provided by dbCAN

PROTEINS="$1"
PREFIX="$2"
DBCAN_HMM="$3"
PARSER="$4"

OUTDIR="${PREFIX}_CAZyme"
mkdir -p "${OUTDIR}"

for FILE in "${PROTEINS}"/*_protein_filtered.faa; do
    BASENAME=$(basename "${FILE}" _protein_filtered.faa)
    SAMPLE_OUT="${OUTDIR}/${BASENAME}"
    mkdir -p "${SAMPLE_OUT}"

    # Run hmmscan
    hmmscan --cpu 8 \
        --domtblout "${SAMPLE_OUT}/${BASENAME}.out.dm" \
        "${DBCAN_HMM}" \
        "${FILE}" > "${SAMPLE_OUT}/${BASENAME}.out"

    # Parse results with dbCAN parser
    bash "${PARSER}" "${SAMPLE_OUT}/${BASENAME}.out.dm" > "${SAMPLE_OUT}/${BASENAME}.out.dm.ps"

    # Apply stringent filter (e-value < 1e-15 and coverage > 0.35)
    awk '$5<1e-15 && $10>0.35' \
        "${SAMPLE_OUT}/${BASENAME}.out.dm.ps" \
        > "${SAMPLE_OUT}/${BASENAME}.out.dm.ps.stringent"
done
