library(tidyverse)
library(dplyr)
library(phyloseq)
library(vegan)
library(pheatmap)
library(RColorBrewer)
set.seed(1234)

# ------------------------------------------
# 1. RPKM matrix (ORF × sample)
# ------------------------------------------
df <- read_table("RPKM_matrix.tsv") %>% 
  column_to_rownames("Gene") %>% 
  as.matrix()

# ------------------------------------------
# 2. CAZy annotation
# ------------------------------------------
cazy <- read_table("HMM_cazy.out.dm.ps.stringent", col_names = FALSE)

df_tax <- df[cazy$X3, ]

cazy2 <- as.data.frame(make.unique(cazy$X3, sep = "."))
rownames(df_tax) <- cazy2[,1]

cazy <- as.data.frame(cbind(cazy = cazy$X1, gene = cazy2[,1]))
rownames(cazy) <- cazy2[,1]

# keep CAZy family level (GHxx, CEyy, etc.)
cazy$cazy <- sub("^([A-Z]+[0-9]+(?:_[0-9]+)?).*", "\\1", cazy$cazy)

# ------------------------------------------
# 3. phyloseq object
# ------------------------------------------
Gene <- otu_table(df_tax, taxa_are_rows = TRUE)
TAX  <- tax_table(as.matrix(cazy))

meta <- read_delim(
  "sample-metadata.txt",
  delim = "\t"
) %>%
  as.data.frame() %>%
  column_to_rownames("ID")

physeq <- phyloseq(Gene, TAX, sample_data(meta))

# ------------------------------------------
# 4. Aggregate to CAZy family
# ------------------------------------------
Gene.cazy <- tax_glom(physeq, taxrank = "cazy")

Gene.cazy.sub <- subset_taxa(
  Gene.cazy,
  grepl("^(CBM|CE|PL|GH)", cazy)
)

comm <- t(as(otu_table(Gene.cazy.sub), "matrix"))

# ======================================================
# 5. CAZy abundance table (for stats & plots)
# ======================================================
cazy_map <- cazy %>%
  dplyr::select(gene, cazy) %>%
  filter(!is.na(cazy))

comm_long <- as.data.frame(comm) %>%
  tibble::rownames_to_column("SampleID") %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "gene",
    values_to = "value"
  ) %>%
  left_join(cazy_map, by = "gene") %>%
  filter(!is.na(cazy))

comm_cazy_long <- comm_long %>%
  group_by(SampleID, cazy) %>%
  summarise(value = sum(value), .groups = "drop")

comm_cazy <- comm_cazy_long %>%
  pivot_wider(
    names_from = cazy,
    values_from = value,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(comm_cazy) <- comm_cazy$SampleID
comm_cazy$SampleID <- NULL

# ------------------------------------------
# 6. Core CAZy (≥50% samples)
# ------------------------------------------
presence_rate <- colMeans(comm_cazy > 0)
core_cazy <- names(presence_rate[presence_rate >= 0.5])

comm_cazy <- comm_cazy[, core_cazy, drop = FALSE]

# ------------------------------------------
# 7. Kruskal–Wallis test
# ------------------------------------------
df_cazy <- as.data.frame(comm_cazy)
df_cazy$Breed <- meta[rownames(df_cazy), "Breed"]

cazy_res <- map_dfr(colnames(comm_cazy), function(cz){
  
  x <- df_cazy %>%
    filter(Breed == "F1") %>%
    pull(.data[[cz]])
  
  y <- df_cazy %>%
    filter(Breed == "JB") %>%
    pull(.data[[cz]])
  
  mean_F1 <- mean(x, na.rm = TRUE)
  mean_JB <- mean(y, na.rm = TRUE)
  
  log2FC <- log2((mean_F1 + 1e-6) / (mean_JB + 1e-6))
  
  wt <- suppressWarnings(
    wilcox.test(x, y, exact = FALSE)
  )
  
  tibble(
    CAZy = cz,
    JB_mean = mean_JB,
    F1_mean = mean_F1,
    log2FC = log2FC,
    W = unname(wt$statistic),
    p = wt$p.value
  )
  
}) %>%
  mutate(FDR = p.adjust(p, method = "fdr")) %>%
  arrange(FDR)

sig_cazy <- cazy_res %>%
  filter(FDR < 0.05)

sig_names <- sig_cazy$CAZy

write.csv(sig_cazy, "Significant_CAZy_statistics.csv", row.names = FALSE)

# ------------------------------------------
# 8. Heatmap (log10 + Z-score)
# ------------------------------------------
mat_sig <- comm_cazy[, sig_names, drop = FALSE]
mat_log <- log10(mat_sig + 1)
mat_z <- scale(mat_log)

dir_df <- df_cazy %>%
  select(all_of(sig_names), Breed) %>%
  pivot_longer(-Breed, names_to = "CAZy", values_to = "value") %>%
  group_by(CAZy, Breed) %>%
  summarise(median = median(value), .groups = "drop") %>%
  pivot_wider(names_from = Breed, values_from = median) %>%
  mutate(Direction = ifelse(F1 > JB, "F1", "JB"))

anno_col <- dir_df %>%
  select(CAZy, Direction) %>%
  as.data.frame()
rownames(anno_col) <- anno_col$CAZy
anno_col$CAZy <- NULL

anno_row <- data.frame(Breed = meta$Breed)
rownames(anno_row) <- rownames(meta)

tiff("Heatmap_CAZy.tiff",
     height = 2600, width = 1800,
     res = 300, compression = "lzw")

pheatmap(
  t(mat_z),
  annotation_row = anno_col,
  annotation_col = anno_row,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-2, 2, length.out = 101),
  cluster_rows = T,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  fontsize = 8
)

dev.off()

# ------------------------------------------
# 9. Bubble plot
# ------------------------------------------
eps <- 1e-6

dir_df <- df_cazy %>%
  select(all_of(sig_names), Breed) %>%
  pivot_longer(-Breed, names_to = "CAZy", values_to = "value") %>%
  group_by(CAZy, Breed) %>%
  summarise(mean_value = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = Breed, values_from = mean_value) %>%
  mutate(
    log2FC = log2((F1 + eps) / (JB + eps)),
    Direction = ifelse(log2FC > 0, "F1", "JB")
  )

vol_df <- dir_df %>%
  left_join(sig_cazy, by = "CAZy")

write.csv(vol_df, "Cazy.csv")

tiff("Bubble_CAZy.tiff",
     height = 3500, width = 1000,
     res = 300, compression = "lzw")

ggplot(vol_df, aes(
  x = Direction,
  y = CAZy,
  color = FDR,
  size = abs(log2FC)
)) +
  guides(
    color = guide_colorbar(order = 1),
    size  = guide_legend(order = 2)
  )+
  geom_point(alpha = 0.8) +
  scale_size_continuous(name = "|log2FC|") +
  scale_color_viridis_c(
    direction = 1,
    name = "FDR"
  ) +
  theme_bw() +
  labs(x = "", y = "") +
  theme(
    axis.text = element_text(size = 10.5, colour = "black")
  )

dev.off()
