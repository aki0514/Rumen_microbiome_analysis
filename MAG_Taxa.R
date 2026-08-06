library(tidyverse)

mag_tax <- read_csv("MAGtaxa.csv")

mag_F1 <- mag_tax %>%
  filter(Breed == "F1")
mag_F1 <- mag_F1 %>%
  mutate(
    Phylum = str_remove(Phylum, "^p__"),
    Genus  = str_remove(Genus,  "^[a-z]__")
  )


phylum_F1 <- mag_F1 %>%
  count(Phylum, name = "MAG_count") %>%
  arrange(desc(MAG_count)) %>%
  mutate(
    Phylum = fct_reorder(Phylum, MAG_count)
  )
tiff("Phylum_MAG.tiff", height = 1600, width = 1600, res = 300, compression = "lzw")
ggplot(phylum_F1, aes(x = Phylum, y = MAG_count, fill = Phylum)) +
  geom_col(width = 0.8) +
  coord_flip() +
  theme_classic() +
  labs(
    x = "Phylum",
    y = "Number of MAGs"
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 12)
  )
dev.off()

genus_F1 <- mag_F1 %>%
  count(Phylum, Genus, name = "MAG_count") %>%
  arrange(desc(MAG_count)) %>%
  slice_head(n = 60) %>%
  mutate(
    Genus_clean = str_remove(Genus, "^[a-z]__"),
    Genus_label = paste0("italic('", Genus_clean, "')"),
    Genus_label = factor(Genus_label, levels = rev(Genus_label))
  )

tiff("Genus_MAG.tiff", height = 2400, width = 1800, res = 300, compression = "lzw")
ggplot(genus_F1, aes(x = Genus_label, y = MAG_count, fill = Phylum)) +
  geom_col(width = 0.8, color = NA) +
  coord_flip() +
  scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  theme_classic() +
  labs(
    x = "Genus",
    y = "Number of MAGs"
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.title = element_blank()
  )
dev.off()

