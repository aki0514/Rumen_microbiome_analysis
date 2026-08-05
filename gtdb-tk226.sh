#!/bin/bash
# Usage:
#   bash gtdb-tk226.sh <GenomeDir> <OutputPrefix> <GTDB_DB>
#
# <GenomeDir>   : directory with dereplicated genomes (e.g., dRep output)
# <OutputPrefix>: prefix for output directory
# <GTDB_DB>     : path to GTDB-Tk reference data (e.g., release226)

GENOMES="$1"
PREFIX="$2"
GTDB_DB="$3"

OUTDIR="${PREFIX}_GTDBtk"
mkdir -p "${OUTDIR}"

# Make sure the DB path is exported
export GTDBTK_DATA_PATH="${GTDB_DB}"

# Run GTDB-Tk
gtdbtk classify_wf \
  --genome_dir "${GENOMES}" \
  --out_dir "${OUTDIR}" \
  --extension fa \
