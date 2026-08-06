library(readr)
library(dplyr)
library(tidyverse)
library(stringr)
library(pheatmap)
library(RColorBrewer)

cazy_mat <- read.csv("cazy_MAG.csv",
                     row.names = 1,
                     check.names = FALSE)

gh_mat <- cazy_mat[grepl("^GH", rownames(cazy_mat)), ]

mag_df <- read_delim("MAG_trait_correlations.tsv",
                     delim = "\t",
                     trim_ws = TRUE)

sig_mag <- mag_df %>%
  filter(FDR < 0.05)
sig_mag <- distinct(sig_mag, MAG, .keep_all = TRUE)

gh_sig <- gh_mat[, sig_mag$MAG]

gh_sig <- gh_sig[rowSums(gh_sig) > 0, ]

gh_sig <- log2(gh_sig + 1)

new_colnames <- paste0(sig_mag$ID, " (", sig_mag$Genus, ")")
new_colnames <- make.unique(new_colnames)
colnames(gh_sig) <- new_colnames

annotation_col <- data.frame(
  Phylum = sig_mag$Phylum,
  Genus  = sig_mag$Genus
)

rownames(annotation_col) <- colnames(gh_sig)

annotation_col$Genus  <- factor(annotation_col$Genus)
annotation_col$Phylum <- factor(annotation_col$Phylum)


phylum_levels <- levels(annotation_col$Phylum)
phylum_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Dark2"))(length(phylum_levels)),
  phylum_levels
)

genus_levels <- levels(annotation_col$Genus)
genus_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set1"))(length(genus_levels)),
  genus_levels
)

ann_colors <- list(
  Phylum = phylum_colors,
  Genus  = genus_colors
)

tiff("GH_Carcass.tiff",
     height = 4600,
     width = 2200,
     res = 300,
     compression = "lzw")

pheatmap(gh_sig,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white","firebrick3"))(50),
         border_color = NA,
         fontsize_col = 6,
         fontsize_row = 6,
         angle_col = 90)

dev.off()
