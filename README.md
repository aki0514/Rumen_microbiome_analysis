# Genome-Resolved Metagenomics Reveals the Taxonomic and Functional Landscape of the Rumen Microbiome in Japanese Black × Holstein Cattle and its Associations with Carcass Traits

This repository contains the shell scripts and R scripts used for the genome-resolved metagenomic analyses described in our manuscript:

> **Genome-Resolved Metagenomics Reveals the Taxonomic and Functional Landscape of the Rumen Microbiome in Japanese Black × Holstein Cattle and its Associations with Carcass Traits**

## Repository structure

```
.
├── shell/
│   ├── trim.sh
│   ├── Assemble.sh
│   ├── BWA_back.sh
│   ├── BWA_reference.sh
│   ├── CheckM2.sh
│   ├── CONCOCT.sh
│   ├── DAS_Tool.sh
│   ├── gtdb-tk226.sh
│   ├── MaxBin2.sh
│   ├── metabat2.sh
│   ├── prodigal.sh
│   ├── cd_hit.sh
│   ├── dRep95.sh
│   ├── dRep99.sh
│   ├── coverm_MAG.sh
│   ├── coverm_protein.sh
│   ├── kofamscan.sh
│   └── HMM_cazy.sh
│
├── R/
│   ├── MAG_diversity.R
│   ├── MAG_abundance.R
│   ├── MAG_Taxa.R
│   ├── MAG_Cor.R
│   ├── KO.R
│   ├── Cazy.R
│   ├── Protein_share.R
│   ├── MAG_cazy.R
│   └── MAG_cazy_carcass.R
```

## Description

### shell/

Shell scripts used for genome-resolved metagenomic analyses, including:

- Gene prediction (Prodigal)
- Protein clustering (CD-HIT)
- Genome dereplication (dRep)
- Genome abundance estimation (CoverM)
- Protein abundance estimation (CoverM)
- Functional annotation (KOfamScan)
- CAZyme annotation (dbCAN/HMMER)

### R/

R scripts used for downstream statistical analyses and figure generation, including:

- Alpha- and beta-diversity analyses
- Differential abundance analyses
- Taxonomic composition
- Correlation analyses with carcass traits
- KEGG pathway analyses
- CAZyme analyses
- Protein repertoire analyses
- Heatmap generation

## Citation

If you use these scripts, please cite:

Sato et al. *Genome-Resolved Metagenomics Reveals the Taxonomic and Functional Landscape of the Rumen Microbiome in Japanese Black × Holstein Cattle and its Associations with Carcass Traits*.
(In preparation)
