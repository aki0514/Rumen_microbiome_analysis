library(tidyverse)
library(pheatmap)
library(RColorBrewer)

cazy_mat <- read.csv("cazy_MAG.csv",
                     row.names = 1,
                     check.names = FALSE)

target_GH <- c("GH5", "GH9", "GH10", "GH11",
               "GH30", "GH43", "GH13", "GH31")

pattern <- paste0("^(", paste(target_GH, collapse="|"), ")(_|\\b)")

gh_mat <- cazy_mat[grepl(pattern, rownames(cazy_mat)), ]

gh_mat <- gh_mat[rowSums(gh_mat) > 0, ]

gh_log <- log1p(gh_mat)

mag_df <- read.csv("Significant_MAGs_with_taxa_andDirection.csv",
                   stringsAsFactors = FALSE)

sig_mag <- mag_df %>%
  filter(FDR < 0.05)

sig_mag$Genus <- sub("^g__", "", sig_mag$Genus)
sig_mag$Genus[sig_mag$Genus == "" | is.na(sig_mag$Genus)] <- "Unknown"

genus_counts <- sig_mag %>%
  count(Genus)

sig_mag <- sig_mag %>%
  left_join(genus_counts, by = "Genus") %>%
  mutate(Genus2 = ifelse(n == 1, "Others", Genus))

sig_mag$Higher <- factor(sig_mag$Higher, levels = c("JB", "F1"))

sig_mag <- sig_mag %>%
  arrange(Higher, Genus2)

gh_sig <- gh_log[, sig_mag$MAG]

new_colnames <- paste0(sig_mag$ID, " (", sig_mag$Genus, ")")
colnames(gh_sig) <- new_colnames

annotation_col <- data.frame(
  Breed = sig_mag$Higher,
  Genus = sig_mag$Genus2
)

rownames(annotation_col) <- colnames(gh_sig)

annotation_col$Genus <- factor(annotation_col$Genus)

genus_levels <- levels(annotation_col$Genus)

genus_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set1"))(length(genus_levels)),
  genus_levels
)

if ("Others" %in% genus_levels) {
  genus_colors["Others"] <- "grey70"
}

ann_colors <- list(
  Breed = c(JB = "#1f78b4", F1 = "#e31a1c"),
  Genus = genus_colors
)

tiff("GH_JBvsF1.tiff", height = 2800, width = 3500, res = 300, compression = "lzw")
pheatmap(gh_sig,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white","firebrick3"))(50),
         border_color = NA,
         fontsize_col = 6,
         fontsize_row = 6,
         angle_col = 90)
dev.off()
