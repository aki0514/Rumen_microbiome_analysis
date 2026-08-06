library(tidyverse)
library(phyloseq)
library(vegan)
set.seed(1234)

# ------------------------------------------
# 1. Read RPKM matrix (ORF × Sample)
# ------------------------------------------
df_rpkm <- read_table("RPKM_matrix.tsv") %>%
  as.data.frame()

# ------------------------------------------
# 2. Read KO annotation (kofamscan)
# ------------------------------------------
ko <- read_tsv(
  "F1_kofamscan.txt",
  col_names = c("Gene", "KO")
) %>%
  filter(!is.na(KO))

# ------------------------------------------
# 3. Join RPKM with KO (safe, ordered)
# ------------------------------------------
df_rpkm_ko <- df_rpkm %>%
  inner_join(ko, by = "Gene")

# ------------------------------------------
# 4. Long format (ORF → KO → Sample)
# ------------------------------------------
df_long <- df_rpkm_ko %>%
  pivot_longer(
    cols = -c(Gene, KO),
    names_to  = "SampleID",
    values_to = "value"
  )

# ------------------------------------------
# 5. Aggregate ORFs at KO level (Sample × KO)
# ------------------------------------------
df_ko_long <- df_long %>%
  group_by(SampleID, KO) %>%
  summarise(value = sum(value), .groups = "drop")

df_ko <- df_ko_long %>%
  pivot_wider(
    names_from  = KO,
    values_from = value,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(df_ko) <- df_ko$SampleID
df_ko$SampleID  <- NULL

# ------------------------------------------
# 7. Metadata
# ------------------------------------------
meta <- read_delim(
  "sample-metadata.txt",
  delim = "\t"
) %>%
  as.data.frame() %>%
  column_to_rownames("ID")

df_ko <- df_ko[rownames(meta), ]

# ------------------------------------------
# 8. Hellinger transformation
# ------------------------------------------
# Hellinger = sqrt(relative abundance)
comm_hel <- decostand(
  df_ko,
  method = "hellinger"
)

# ------------------------------------------
# 9. NMDS (Euclidean is correct for Hellinger)
# ------------------------------------------
ord.nmds <- metaMDS(
  comm_hel,
  distance = "euclidean",
  k = 2,
  trymax = 200,
  autotransform = FALSE
)

ord.nmds$stress   # 必ず確認

# ------------------------------------------
# 10. NMDS plot
# ------------------------------------------
adonis_res <- adonis2(
  comm_hel ~ Breed,
  data = meta,
  method = "euclidean",
  permutations = 999
)

r2 <- adonis_res$R2[1]
p  <- adonis_res$`Pr(>F)`[1]

lab <- paste0(
  "R² = ", sprintf("%.3f", r2),
  ", P ", ifelse(p < 0.05, "< 0.05", sprintf("%.3f", p))
)

df_nmds <- as.data.frame(ord.nmds$points)
df_nmds$SampleID <- rownames(df_nmds)
df_nmds$Breed <- meta[df_nmds$SampleID, "Breed"]

tiff("NMDS_KO_Hellinger.tiff",
     height = 1300, width = 1400,
     res = 300, compression = "lzw")

ggplot(df_nmds, aes(x = MDS1, y = MDS2, color = Breed)) +
  geom_point(size = 3) +
  stat_ellipse(type = "t", linetype = 2, linewidth = 0.5)+
  annotate(
    "text",
    x = Inf, y = 0.25,
    label = lab,
    hjust = 1.1, vjust = 0.5,
    size = 3.5
  ) +
  scale_color_manual(values = c(JB="#1f78b4", F1="#e31a1c")) +
  theme_bw() +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    title = ""
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10, color = "black")
  )

dev.off()

# ============================
# ============================

dist_bc <- vegdist(df_ko, method = "bray")
meta <- read_tsv("sample-metadata.txt")

# ============================
# 11) PERMANOVA（BC）
# ============================

adonis_bc <- adonis2(
  dist_bc ~ Breed,
  data = meta,
  permutations = 999
)

print(adonis_bc)

# ============================
# 12) NMDS（BC）
# ============================

nmds_bc <- metaMDS(
  dist_bc,
  k = 2,
  trymax = 100
)

nmds_bc$stress

# ============================
# ============================

nmds_bc_points <- as.data.frame(nmds_bc$points) %>%
  rownames_to_column("ID") %>%
  left_join(meta, by = "ID")

r2_bc <- adonis_bc$R2[1]
p_bc  <- adonis_bc$`Pr(>F)`[1]

lab_bc <- paste0(
  "R² = ", sprintf("%.3f", r2_bc),
  ", P ", ifelse(p_bc < 0.05, "< 0.05", sprintf("%.3f", p_bc))
)

# ============================
# ============================

tiff(
  "NMDS_KO_BrayCurtis_Breed.tiff",
  height = 1300,
  width  = 1500,
  res    = 300,
  compression = "lzw"
)

ggplot(nmds_bc_points, aes(x = MDS1, y = MDS2, color = Breed)) +
  geom_point(size = 3) +
  stat_ellipse(type = "t", linetype = 2, linewidth = 0.5) +
  annotate(
    "text",
    x = Inf, y = 0.26,
    label = lab_bc,
    hjust = 1.1, vjust = 0.5,
    size = 3.5
  ) +
  scale_color_manual(values = c(JB = "#1f78b4", F1 = "#e31a1c")) +
  theme_bw() +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    title = ""
  ) +
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    panel.grid = element_blank()
  )

dev.off()

# ------------------------------------------
# 11. PERMANOVA (same transformed data)
# ------------------------------------------


ko2path <- read_tsv(
  "ko2pathway.tsv",
  col_names = c("KO", "Pathway")
) %>%
  filter(str_detect(Pathway, "^ko"))

df_ko_long2 <- df_ko %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "KO",
    values_to = "value"
  )

df_path_long <- df_ko_long2 %>%
  left_join(ko2path, by = "KO") %>%
  filter(!is.na(Pathway))

df_path <- df_path_long %>%
  group_by(SampleID, Pathway) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  pivot_wider(
    names_from = Pathway,
    values_from = value,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(df_path) <- df_path$SampleID
df_path$SampleID <- NULL

df_path_stat <- df_path %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(
    meta %>% rownames_to_column("SampleID"),
    by = "SampleID"
  )

df_path_stat <- df_path_stat[, 1:(ncol(df_path_stat)-2)]

path_long <- df_path_stat %>%
  pivot_longer(
    cols = -c(SampleID, Breed),
    names_to = "Pathway",
    values_to = "value"
  )

eps <- 1e-6

path_mean <- path_long %>%
  group_by(Pathway, Breed) %>%
  summarise(mean_value = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = Breed, values_from = mean_value) %>%
  mutate(
    log2FC = log2((F1 + eps) / (JB + eps)),
    Direction = ifelse(log2FC > 0, "F1", "JB")
  )

path_test <- path_long %>%
  group_by(Pathway) %>%
  summarise(
    W = unname(wilcox.test(value ~ Breed, exact = FALSE)$statistic),
    p = wilcox.test(value ~ Breed, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(p, method = "fdr"))

vol_path <- path_mean %>%
  left_join(path_test, by = "Pathway")

sig_path <- vol_path %>%
  filter(FDR < 0.05) %>%
  select(
    Pathway,
    JB,
    F1,
    log2FC,
    W,
    p,
    FDR
  ) %>%
  arrange(FDR)

write.csv(sig_path,
          "Significant_Pathways.csv",
          row.names = FALSE)

tiff("Bubble_Pathway.tiff",
     height = 2200, width = 1100,
     res = 300, compression = "lzw")

ggplot(vol_path %>% filter(FDR < 0.05),
       aes(
         x = Direction,
         y = Pathway,
         size = abs(log2FC),
         color = FDR
       )) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_c(
    direction = 1,
    name = "FDR"
  ) +
  labs(
    x = "",
    y = "",
    title = ""
  )+
  scale_size_continuous(name = "|log2FC|") +
  theme_bw() +
  theme(axis.text = element_text(color = "black"))

dev.off()

