# Relationships between rumen microbiome and carcass traits in F1 cattle

This repository contains the shell scripts and R scripts used for the genome-resolved metagenomic analyses described in our manuscript:

> **Relationships between rumen microbiome and carcass traits in F1 cattle**

## Repository structure

```
.
├── shell/
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
│   ├── JBvsF1.R
│   ├── Taxa.R
│   ├── Cor.R
│   ├── KO.R
│   ├── Cazy.R
│   ├── Protein_share.R
│   ├── GH_JBvsF1.R
│   └── cazy_carcass.R
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

Sato et al. *Relationships between rumen microbiome and carcass traits in F1 cattle*.
(In preparation)
