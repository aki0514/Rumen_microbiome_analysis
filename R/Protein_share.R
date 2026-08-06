
library(tidyverse)
library(eulerr)
library(tidytext)
library(ggplot2)
############################################
## 1. Read CD-HIT cluster file
############################################
clstr <- readLines("Merge_protein_nr.faa.clstr")

cluster_df <- tibble(line = clstr) %>%
  mutate(
    cluster = if_else(
      str_detect(line, "^>Cluster"),
      str_remove(line, "^>Cluster "),
      NA_character_
    )
  ) %>%
  fill(cluster)

############################################
## 2. Extract cluster members
############################################
members <- cluster_df %>%
  filter(!str_detect(line, "^>Cluster")) %>%
  mutate(
    cdhit_id = str_extract(line, ">[^|]+") %>% str_remove(">"),
    is_rep   = str_detect(line, "\\*"),
    Breed = case_when(
      str_detect(cdhit_id, "^F1_") ~ "F1",
      str_detect(cdhit_id, "^JB_") ~ "JB",
      TRUE ~ "Other"
    )
  ) %>%
  select(cluster, cdhit_id, Breed, is_rep)

############################################
## 3. Define cluster category
############################################
cluster_class <- members %>%
  group_by(cluster) %>%
  summarise(
    has_F1 = any(Breed == "F1"),
    has_JB = any(Breed == "JB"),
    .groups = "drop"
  ) %>%
  mutate(
    Category = case_when(
      has_F1 & has_JB ~ "Shared",
      has_F1          ~ "F1_only",
      has_JB          ~ "JB_only"
    )
  )

table(cluster_class$Category)

############################################
## 4. Euler diagram (ORF clusters)
############################################
cluster_counts <- table(cluster_class$Category)

fit <- euler(c( "F1" = 1188519 + 1111463, "JB" = 1616364 + 1111463, "F1&JB" = 1111463 ))

tiff("Venn_ORF_clusters.tiff",
     height = 2600, width = 1800,
     res = 300, compression = "lzw")

plot(
  fit,
  fills = c("#E69F00", "#56B4E9"),
  edges = TRUE,
  labels = list(font = 2, cex = 1.2),
  quantities = list(font = 2, cex = 1.1)
)

dev.off()

############################################
## 5. Representative ORF per cluster
############################################
rep_orf <- members %>%
  filter(is_rep) %>%
  select(cluster, cdhit_id, Breed)

############################################
## 6. Map CD-HIT ID to original ORF ID
############################################
id_map <- readLines("ID.txt") %>%
  tibble(line = .) %>%
  mutate(
    line     = str_remove(line, "^>"),
    cdhit_id = str_extract(line, "^[^|]+"),
    kofam_id = str_extract(line, "(?<=\\|).+$")
  ) %>%
  select(cdhit_id, kofam_id)

############################################
## 7. Read KOfamScan results
############################################
ko <- read_tsv(
  "F1_kofamscan.txt",
  col_names = c("kofam_id", "KO")
)

############################################
## 8. Assign KO to representative ORFs
############################################
rep_orf_ko <- rep_orf %>%
  left_join(id_map, by = "cdhit_id") %>%
  left_join(ko, by = "kofam_id")

############################################
## 9. Cluster × KO × Category table
############################################
cluster_ko_all <- rep_orf_ko %>%
  left_join(cluster_class, by = "cluster") %>%
  filter(!is.na(KO)) %>%
  distinct(cluster, KO, Category)

############################################
## 10. Count KO by category (cluster-based)
############################################
ko_by_category <- cluster_ko_all %>%
  count(Category, KO, sort = TRUE)

############################################
## 11. Top 30 KO per category
############################################
top_ko <- ko_by_category %>%
  group_by(Category) %>%
  slice_max(n, n = 30) %>%
  ungroup()

tiff("Top30_KO_by_category.tiff",
     height = 1900, width = 1900,
     res = 300, compression = "lzw")

ggplot(top_ko,
       aes(
         x = reorder_within(KO, n, Category),
         y = n,
         fill = Category
       )) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~Category, scales = "free_y") +
  coord_flip() +
  scale_x_reordered() +
  labs(
    x = "",
    y = "Count",
    title = ""
  )+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black")
  )

dev.off()
