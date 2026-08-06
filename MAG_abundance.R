library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(pheatmap)
library(RColorBrewer)

# ----------------------------
# ----------------------------
detect_thresh  <- 0.0   # 検出下限（例 0, 0.05, 0.1）
prev_cut  <- 0.5   # 対称コア: 両群で ≥ 50%
alpha_fdr      <- 0.05  # 有意閾値

rpkm <- read_tsv("RPKM_matrix.tsv") %>% column_to_rownames("MAG")  
meta <- read_tsv("sample-metadata.txt")                         

rpkm_t <- as.data.frame(t(rpkm)) %>% rownames_to_column("ID")
df <- left_join(meta, rpkm_t, by = "ID")

df$Breed <- factor(df$Breed, levels = c("JB","F1"))
JB_ids <- df %>% filter(Breed=="JB") %>% pull(ID)
F1_ids <- df %>% filter(Breed=="F1") %>% pull(ID)

pres <- (rpkm > detect_thresh) * 1
prev_df <- tibble(
  MAG     = rownames(rpkm),
  prev_JB = rowMeans(pres[, intersect(JB_ids, colnames(pres)), drop=FALSE], na.rm=TRUE),
  prev_F1 = rowMeans(pres[, intersect(F1_ids, colnames(pres)), drop=FALSE], na.rm=TRUE)
)

core_asym <- prev_df %>%
  filter(prev_JB >= prev_cut | prev_F1 >= prev_cut) %>%
  pull(MAG)

rpkm_core <- rpkm[core_asym, , drop=FALSE]
rpkm_core_t <- as.data.frame(t(rpkm_core)) %>%
  rownames_to_column("ID")

df_core <- left_join(meta, rpkm_core_t, by="ID") %>%
  mutate(Breed = factor(Breed, levels=c("JB","F1")))

wilcox_res <- map_dfr(core_asym, function(mag){
  
  x <- df_core %>% filter(Breed=="F1") %>% pull(.data[[mag]])
  y <- df_core %>% filter(Breed=="JB") %>% pull(.data[[mag]])
  
  mean_F1 <- mean(x, na.rm=TRUE)
  mean_JB <- mean(y, na.rm=TRUE)
  
  # log2 Fold Change (F1 / JB)
  log2fc <- log2((mean_F1 + 1e-6) / (mean_JB + 1e-6))
  
  if(sum(!is.na(x)) < 2 || sum(!is.na(y)) < 2){
    
    tibble(
      MAG = mag,
      JB_mean = mean_JB,
      F1_mean = mean_F1,
      W = NA_real_,
      pval = NA_real_,
      log2FC = log2fc
    )
    
  } else {
    
    wt <- suppressWarnings(
      wilcox.test(x, y, exact = FALSE)
    )
    
    tibble(
      MAG = mag,
      JB_mean = mean_JB,
      F1_mean = mean_F1,
      W = unname(wt$statistic),
      pval = wt$p.value,
      log2FC = log2fc
    )
    
  }
  
}) %>%
  mutate(FDR = p.adjust(pval, method = "fdr"))

res2 <- wilcox_res %>%
  left_join(prev_df, by="MAG") %>%
  select(
    MAG,
    JB_mean,
    F1_mean,
    prev_JB,
    prev_F1,
    log2FC,
    W,
    pval,
    FDR
  ) %>%
  arrange(FDR)

write_tsv(res2, "CoreMAG_diff_JB_vs_F1.tsv")

# ----------------------------
# ----------------------------
sig_top <- res2 %>%
  filter(!is.na(FDR), FDR < alpha_fdr) %>%
  arrange(FDR) %>%
  pull(MAG)

mat_sig <- rpkm_core[intersect(sig_top, rownames(rpkm_core)), , drop = FALSE]
sample_order <- c(JB_ids, F1_ids) %>% intersect(colnames(mat_sig))
mat_sig <- mat_sig[, sample_order, drop = FALSE]

# --- log + Z-score ---
mat_scaled <- t(scale(t(log1p(as.matrix(mat_sig)))))

# ----------------------------
# column annotation（Breed）
# ----------------------------
ann_col <- meta %>%
  select(ID, Breed) %>%
  column_to_rownames("ID")

# ----------------------------
# row annotation（direction + Phylum）
# ----------------------------

JB_cols <- intersect(JB_ids, colnames(mat_sig))
F1_cols <- intersect(F1_ids, colnames(mat_sig))
JB_mean <- rowMeans(mat_sig[, JB_cols, drop = FALSE], na.rm = TRUE)
F1_mean <- rowMeans(mat_sig[, F1_cols, drop = FALSE], na.rm = TRUE)

Direction <- ifelse(JB_mean > F1_mean, "JB", "F1")

taxa <- read_csv("MAGtaxa.csv")

row_taxa <- taxa %>%
  select(MAG, Phylum) %>%
  column_to_rownames("MAG")

row_anno <- data.frame(
  Direction = Direction,
  Phylum    = row_taxa[rownames(mat_sig), "Phylum"],
  row.names = rownames(mat_sig)
)

row_anno$Phylum <- gsub("^p__", "", row_anno$Phylum)
# ----------------------------
# annotation colors
# ----------------------------

phylum_levels <- unique(na.omit(row_anno$Phylum))
phylum_cols <- setNames(
  colorRampPalette(brewer.pal(8, "Set1"))(length(phylum_levels)),
  phylum_levels
)

ann_colors <- list(
  Breed     = c(JB = "#1f78b4", F1 = "#e31a1c"),
  Direction = c(JB = "#1f78b4", F1 = "#e31a1c"),
  Phylum    = phylum_cols
)

gaps_col <- length(JB_ids)

row_anno$Direction <- factor(row_anno$Direction, levels = c("JB", "F1"))

row_order <- order(row_anno$Direction)

mat_scaled <- mat_scaled[row_order, , drop = FALSE]
row_anno   <- row_anno[row_order, , drop = FALSE]

gaps_row <- sum(row_anno$Direction == "JB")

# ----------------------------
# ----------------------------

bk <- seq(-2, 2, length = 101)
hm_cols <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

tiff("Heatmap_CoreMAGs_Breed.tiff",
     width = 2000, height = 2000, res = 300, compression = "lzw")

pheatmap(
  mat_scaled,
  annotation_col    = ann_col,
  annotation_row    = row_anno,
  annotation_colors = ann_colors,
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,   
  gaps_col          = length(JB_ids),
  gaps_row          = gaps_row,
  color             = hm_cols,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  breaks            = bk
)


dev.off()


sig_res <- res2 %>%
  filter(!is.na(FDR), FDR < 0.05) %>%
  arrange(FDR)

taxa <- read_csv("MAGtaxa.csv")

JB_ids <- meta %>% filter(Breed=="JB") %>% pull(ID)
F1_ids <- meta %>% filter(Breed=="F1") %>% pull(ID)

grp_means <- tibble(
  MAG = rownames(rpkm),
  JB  = rowMeans(rpkm[, intersect(JB_ids, colnames(rpkm)), drop=FALSE], na.rm=TRUE),
  F1  = rowMeans(rpkm[, intersect(F1_ids, colnames(rpkm)), drop=FALSE], na.rm=TRUE)
)

sig_res2 <- sig_res %>%
  left_join(grp_means, by="MAG") %>%
  mutate(
    Higher = case_when(
      JB > F1 ~ "JB",
      F1 > JB ~ "F1",
      TRUE    ~ "Equal"
    )
  )

sig_res_taxa <- sig_res2 %>%
  left_join(taxa, by = "MAG")

write.csv(sig_res_taxa, "Significant_MAGs_with_taxa_andDirection.csv")


df_phylum_dir <- row_anno %>%
  as.data.frame() %>%
  rownames_to_column("MAG") %>%
  mutate(
    Phylum = gsub("^p__", "", Phylum)  # 念のため
  )

count_df <- df_phylum_dir %>%
  count(Phylum, Direction, name = "n_MAG")

phylum_order <- count_df %>%
  group_by(Phylum) %>%
  summarise(total = sum(n_MAG)) %>%
  arrange(desc(total)) %>%
  pull(Phylum)

count_df$Phylum <- factor(
  count_df$Phylum,
  levels = rev(phylum_order)
)
count_df$Direction <- factor(count_df$Direction, levels = c("JB","F1"))

p_bar <- ggplot(count_df,
                aes(x = Phylum, y = n_MAG, fill = Direction)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(
    values = c(JB = "#1f78b4", F1 = "#e31a1c")
  ) +
  theme_bw(base_size = 11) +
  labs(
    x = "",
    y = "Number of significant MAGs",
    fill = ""
  ) +
  theme(axis.text = element_text(size = 10, colour = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

tiff("Barplot_Phylum_direction_MAGcount.tiff",
     width = 1400, height = 1400, res = 300, compression = "lzw")
print(p_bar)
dev.off()


library(forcats)

sig_genus <- sig_res_taxa %>%
  mutate(
    Genus = gsub("^g__", "", Genus),
    Genus = ifelse(is.na(Genus) | Genus == "",
                   "Unclassified", Genus)
  ) %>%
  count(Genus, Higher, name = "n_MAG") 

top_genus <- sig_genus %>%
  group_by(Genus) %>%
  summarise(total = sum(n_MAG)) %>%
  arrange(desc(total)) %>%
  slice_head(n = 100) %>%
  pull(Genus)

plot_df <- sig_genus %>%
  filter(Genus %in% top_genus) %>%
  mutate(
    Genus = fct_reorder(Genus, n_MAG, .fun = sum, .desc = FALSE)
  )

tiff("Barplot_Genus_direction_MAGcount.tiff",
     width = 1400, height = 3000, res = 300, compression = "lzw")

ggplot(plot_df, aes(x = Genus, y = n_MAG, fill = Higher)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c(JB="#1f78b4", F1="#e31a1c")) +
  scale_y_continuous(
    breaks = function(x) seq(0, ceiling(max(x)), by = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_bw(base_size = 13) +
  labs(
    x = "",
    y = "Number of MAGs",
    fill = ""
  ) +
  theme(axis.text = element_text(size = 10, colour = "black"),
    axis.text.y = element_text(face = "italic"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

dev.off()
