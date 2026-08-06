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
# 6. Metadata
# ------------------------------------------
meta <- read_delim(
  "sample-metadata.txt",
  delim = "\t"
) %>%
  as.data.frame() %>%
  column_to_rownames("ID")

df_ko <- df_ko[rownames(meta), ]

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

