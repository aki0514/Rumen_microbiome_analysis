library(vegan)
library(tidyverse)
library(ape)
library(ggplot2)
set.seed(123)

rpkm_mat <- read_tsv("RPKM_matrix.tsv") %>%
  column_to_rownames("MAG") %>%
  as.matrix()

meta <- read_tsv("sample-metadata.txt")
meta <- meta %>%
  filter(ID %in% colnames(rpkm_mat)) %>%
  arrange(match(ID, colnames(rpkm_mat)))

rpkm_hel <- decostand(t(rpkm_mat), method = "hellinger")
dist_euc <- dist(rpkm_hel, method = "euclidean")
adonis_result <- adonis2(
  dist_euc ~ Breed,
  data = meta,
  permutations = 999
)

print(adonis_result)

disp_euc <- betadisper(dist_euc, meta$Breed)

anova(disp_euc)

permutest(disp_euc, permutations = 999)

nmds <- metaMDS(
  dist_euc,
  k = 2,
  trymax = 100,
  autotransform = FALSE
)

nmds_points <- as.data.frame(nmds$points) %>%
  rownames_to_column("ID") %>%
  left_join(meta, by = "ID")

r2 <- adonis_result$R2[1]
p  <- adonis_result$`Pr(>F)`[1]
lab <- paste0(
  "R² = ", sprintf("%.3f", r2),
  ", P = ", ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
)

tiff(
  "NMDS_MAG_Hellinger_Euclidean_Breed.tiff",
  height = 1300,
  width  = 1500,
  res    = 300,
  compression = "lzw"
)

ggplot(nmds_points, aes(x = MDS1, y = MDS2, color = Breed)) +
  geom_point(size = 3) +
  stat_ellipse(type = "t", linetype = 2, linewidth = 0.5) +
  annotate(
    "text",
    x = Inf, y = 0.35,
    label = lab,
    hjust = 1.1, vjust = 1.5,
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
    axis.title = element_text(size = 11),    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 1),
    legend.background = element_rect(fill = NA, colour = "black"),
    legend.key = element_rect(fill = "white", color = NA)
  )

dev.off()


dist_bc <- vegdist(t(rpkm_mat), method = "bray")
adonis_bc <- adonis2(
  dist_bc ~ Breed,
  data = meta,
  permutations = 999
)

print(adonis_bc)

disp_bc <- betadisper(dist_bc, meta$Breed)

anova(disp_bc)

permutest(disp_bc, permutations = 999)

nmds_bc <- metaMDS(
  dist_bc,
  k = 2,
  trymax = 100
)

nmds_bc_points <- as.data.frame(nmds_bc$points) %>%
  rownames_to_column("ID") %>%
  left_join(meta, by = "ID")

r2_bc <- adonis_bc$R2[1]
p_bc  <- adonis_bc$`Pr(>F)`[1]

lab_bc <- paste0(
  "R² = ", sprintf("%.3f", r2_bc),
  ", P = ", ifelse(p_bc < 0.001, "< 0.001", sprintf("%.3f", p_bc))
)

tiff(
  "NMDS_MAG_BrayCurtis_Breed.tiff",
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
    x = Inf, y = 0.8,
    label = lab_bc,
    hjust = 1.1, vjust = 1.5,
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
    axis.title = element_text(size = 11),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 1),
    legend.background = element_rect(fill = NA, colour = "black"),
    legend.key = element_rect(fill = "white", color = NA)
  )

dev.off()

observed <- colSums(rpkm_mat > 0)

alpha_df <- tibble(
  ID = names(observed),
  Observed = observed
) %>%
  left_join(meta, by = "ID")

shannon <- diversity(t(rpkm_mat), index = "shannon")

alpha_df <- alpha_df %>%
  mutate(Shannon = shannon)

print(wilcox.test(Observed ~ Breed, data = alpha_df))
print(wilcox.test(shannon ~ Breed, data = alpha_df))

tiff("Observed_Breed.tiff", height = 1300, width = 1500, res = 300, compression = "lzw")
ggplot(alpha_df, aes(x = Breed, y = Observed, fill = Breed)) +
  geom_violin(trim = TRUE, alpha=0.5, linewidth=0.2) +
  geom_boxplot(width = 0.25, outlier.shape = NA, alpha=0.7) +
  geom_jitter(width = 0.12, size=1.5, alpha=0.7) +
  scale_fill_manual(values = c(JB="#1f78b4", F1="#e31a1c")) +
  theme_bw() +
  labs(y = "Observed MAGs", x = NULL,
       title = "")+
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 1),
    legend.background = element_rect(fill = NA, colour = "black"),
    legend.key = element_rect(fill = "white", color = NA)
  )
dev.off()

tiff("Shannon_Breed.tiff", height = 1300, width = 1500, res = 300, compression = "lzw")
ggplot(alpha_df, aes(x = Breed, y = shannon, fill = Breed)) +
  geom_violin(trim = TRUE, alpha=0.5, linewidth=0.2) +
  geom_boxplot(width = 0.25, outlier.shape = NA, alpha=0.7) +
  geom_jitter(width = 0.12, size=1.5, alpha=0.7) +
  scale_fill_manual(values = c(JB="#1f78b4", F1="#e31a1c")) +
  theme_bw() +
  labs(y = "Shannon index", x = NULL,
       title = "")+
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 1),
    legend.background = element_rect(fill = NA, colour = "black"),
    legend.key = element_rect(fill = "white", color = NA)
  )
dev.off()
