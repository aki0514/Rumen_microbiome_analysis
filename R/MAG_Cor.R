library(tidyverse)
library(purrr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(tidytext)
`%||%` <- function(a, b) if (!is.null(a)) a else b

rpkm <- read_tsv("RPKM_matrix.tsv", show_col_types = FALSE) %>%
  column_to_rownames("MAG") %>% 
  as.matrix()

meta <- read_tsv("sample-metadata.txt", show_col_types = FALSE) %>%
  filter(Breed == "F1")
colnames(meta)[c(3,5)] <- c("Carcass weight", "Ribeye area")

vfa <- read_csv("VFA.csv", show_col_types = FALSE) %>%
  mutate(across(Acetate:Valerate, ~ .x / Total * 100))

tax <- read_csv("MAGtaxa.csv", show_col_types = FALSE) %>%
  transmute(
    MAG,
    ID,
    Phylum = str_remove(Phylum, "^p__"),
    Genus = str_remove(Genus, "^g__")
  ) %>%
  mutate(
    Genus = if_else(is.na(Genus) | Genus == "", "Unclassified", Genus),
    Label = paste0(ID, " (", Genus, ")")
  )

mag2label <- setNames(tax$Label, tax$MAG)

df <- meta %>%
  left_join(
    as.data.frame(t(rpkm)) %>% rownames_to_column("Sample"),
    by = c("ID" = "Sample")
  )

mag_cols <- setdiff(colnames(df), colnames(meta))
prev <- sapply(mag_cols, \(m) mean(df[[m]] > 0, na.rm = TRUE))
keep_mags <- names(prev)[prev >= 0.5]

calc_corr <- function(mag, trait){
  x <- log1p(df[[mag]])
  y <- df[[trait]]
  ok <- complete.cases(x, y)
  if (sum(ok) < 5) return(NULL)
  
  ct <- cor.test(x[ok], y[ok], method = "spearman")
  tibble(
    MAG = mag,
    Trait = trait,
    Cor = unname(ct$estimate),
    Pval = ct$p.value
  )
}

traits <- c("Marbling", "Carcass weight", "Ribeye area", "Backfat", "Rib_thickness")

res <- map_dfr(keep_mags, \(m){
  map_dfr(traits, ~calc_corr(m, .x))
}) %>%
  group_by(Trait) %>%
  mutate(FDR = p.adjust(Pval, "fdr")) %>%
  ungroup() %>%
  left_join(enframe(prev, name="MAG", value="Prevalence"), by="MAG") %>%
  left_join(tax, by="MAG")

write_tsv(res, "MAG_trait_correlations.tsv")

sig <- res %>% filter(FDR < 0.05)
sig <- sig %>%
  arrange(desc(Cor), FDR)
traits_sig <- c("Carcass weight", "Ribeye area")

for (tr in traits_sig) {
  
  sig_sub <- sig %>% filter(Trait == tr)
  
  plots <- map(seq_len(nrow(sig_sub)), function(i){
    mag <- sig_sub$MAG[i]
    
    ggplot(df, aes(x = .data[[tr]], y = log1p(.data[[mag]]))) +
      geom_point(size = 2, alpha = 0.8) +
      geom_smooth(method = "loess", se = FALSE, color = "black") +
      theme_bw() +
      labs(
        title = mag2label[[mag]] %||% mag,
        subtitle = paste0(
          "ρ = ", round(sig_sub$Cor[i], 2),
          ", FDR = ", signif(sig_sub$FDR[i], 2)
        ),
        x = tr,
        y = "MAG abundance (log-transformed)"
      ) +
      theme(
        axis.text = element_text(size = 10, colour = "black"),
        plot.title = element_text(face = "bold", size = 11),
        panel.grid = element_blank()
      )
  })
  
  tiff(
    paste0("Correlation_", tr, ".tiff"),
    height = 5500, width = 5500,
    res = 300, compression = "lzw"
  )
  
  print(wrap_plots(plots, ncol = 6))
  
  dev.off()
}

genus_count <- sig %>%
  filter(Trait %in% c("Weight", "Rib_eye")) %>%
  mutate(
    Direction = ifelse(Cor > 0, "Positive", "Negative")
  ) %>%
  count(Genus, Phylum, Trait, Direction, name = "n") %>%
  mutate(
    n_signed = ifelse(Direction == "Positive", n, -n)
  )

genus_summary <- genus_count %>%
  group_by(Genus) %>%
  summarise(
    total = sum(abs(n_signed)),
    
    has_weight = any(Trait == "Weight"),
    has_ribeye = any(Trait == "Rib_eye"),
    
    has_positive = any(Direction == "Positive"),
    has_negative = any(Direction == "Negative"),
    
    .groups = "drop"
  ) %>%
  mutate(
    trait_class = case_when(
      has_weight & has_ribeye ~ "Both",
      has_weight ~ "Weight",
      has_ribeye ~ "Rib_eye"
    ),
    
    corr_class = case_when(
      has_positive & has_negative ~ 1,
      has_positive ~ 2,
      has_negative ~ 3
    )
  )

genus_order <- genus_summary %>%
  mutate(
    trait_rank = case_when(
      trait_class == "Both" ~ 1,
      trait_class == "Weight" ~ 2,
      trait_class == "Rib_eye" ~ 3
    )
  ) %>%
  arrange(
    trait_rank,      
    corr_class,      
    desc(total)      
  ) %>%
  pull(Genus)

genus_count$Genus <- factor(genus_count$Genus, levels = rev(genus_order))

genus_count <- genus_count %>%
  group_by(Trait, Genus) %>%
  summarise(
    total_trait = sum(abs(n_signed)),
    sign = ifelse(sum(n_signed) > 0, 1, -1),   
    .groups = "drop"
  )

genus_count <- sig %>%
  filter(Trait %in% c("Weight", "Rib_eye")) %>%
  mutate(
    Direction = ifelse(Cor > 0, "Positive", "Negative")
  ) %>%
  count(Genus, Phylum, Trait, Direction, name = "n") %>%
  mutate(
    n_signed = ifelse(Direction == "Positive", n, -n)
  ) %>%
  left_join(genus_count, by = c("Genus", "Trait"))

genus_count <- genus_count %>%
  mutate(
    order_score = sign * total_trait,
    Genus_ordered = reorder_within(Genus, order_score, Trait)
  )

phylum_levels <- unique(genus_count$Phylum)

phylum_cols <- setNames(
  RColorBrewer::brewer.pal(max(3, length(phylum_levels)), "Set2")[seq_along(phylum_levels)],
  phylum_levels
)

tiff("Genus_enrichment_Phylum_PosNeg_Combined.tiff",
     1500, 2500, res = 300, compression = "lzw")
genus_count$Trait <- recode(genus_count$Trait,
                            "Rib_eye" = "Ribeye area",
                            "Weight" = "Carcass weight")

ggplot(genus_count, aes(x = Genus_ordered, y = n_signed, fill = Phylum)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(values = phylum_cols) +
  
  facet_wrap(~ Trait, ncol = 1, scales = "free_y") +
  
  scale_x_reordered() +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    axis.text.y = element_text(face = "italic"),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid = element_blank()
  ) +
  labs(
    x = "",
    y = "Number of MAGs",
    fill = "Phylum"
  )

dev.off()

ord_df <- res %>%
  filter(Trait %in% c("Weight", "Rib_eye"), FDR < 0.05) %>%
  group_by(MAG) %>%
  summarise(
    Cor_weight = Cor[Trait == "Weight"] %>% first() %||% NA,
    Cor_ribeye = Cor[Trait == "Rib_eye"] %>% first() %||% NA,
    
    has_weight = any(Trait == "Weight"),
    has_ribeye = any(Trait == "Rib_eye"),
    
    .groups = "drop"
  ) %>%
  mutate(
    Cor_main = ifelse(!is.na(Cor_weight), Cor_weight, Cor_ribeye),
    
    Correlation = ifelse(Cor_main > 0, "Positive", "Negative"),
    
    Traits = case_when(
      has_weight & has_ribeye ~ "Both",
      has_weight ~ "Weight",
      has_ribeye ~ "Rib_eye"
    )
  ) %>%
  arrange(
    Correlation,
    desc(abs(Cor_main))
  )

ord <- ord_df$MAG
row_anno <- data.frame(
  Correlation = ord_df$Correlation,
  Traits  = ord_df$Traits,
  row.names   = ord
)

ann_cols <- list(
  Correlation = c(
    Positive = "#e31a1c",   # 赤
    Negative = "#1f78b4"    # 青
  ),
  Traits = c(
    'Carcass weight'  = "#1f78b4",    
    'Ribeye area' = "#33a02c",   
    Both    = "#ff7f00"     
  )
)
# ---------------------------------------------------------
# ---------------------------------------------------------
rpkm_t <- as.data.frame(t(rpkm[ord, , drop = FALSE])) %>%
  rownames_to_column("ID")

df_vfa <- meta %>%
  left_join(rpkm_t, by = "ID") %>%
  left_join(vfa, by = "ID")

# ---------------------------------------------------------
# ---------------------------------------------------------
vfa_targets <- c("Acetate", "Propionate", "Butyrate")

calc_cor_vfa <- function(mag, v){
  ct <- cor.test(
    log1p(df_vfa[[mag]]),
    df_vfa[[v]],
    method = "spearman",
    exact = FALSE
  )
  tibble(
    MAG  = mag,
    VFA  = v,
    rho  = unname(ct$estimate),
    pval = ct$p.value
  )
}

cor_vfa <- map_dfr(ord, \(m){
  map_dfr(vfa_targets, ~calc_cor_vfa(m, .x))
}) %>%
  group_by(VFA) %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  ungroup()

# ---------------------------------------------------------
# ---------------------------------------------------------
mat_rho <- cor_vfa %>%
  select(MAG, VFA, rho) %>%
  pivot_wider(names_from = VFA, values_from = rho) %>%
  column_to_rownames("MAG") %>%
  as.matrix()

mat_fdr <- cor_vfa %>%
  select(MAG, VFA, FDR) %>%
  pivot_wider(names_from = VFA, values_from = FDR) %>%
  column_to_rownames("MAG") %>%
  as.matrix()

mat_rho <- mat_rho[ord, , drop = FALSE]
mat_fdr <- mat_fdr[ord, , drop = FALSE]

# ---------------------------------------------------------
# ---------------------------------------------------------
sig_label <- ifelse(mat_fdr < 0.05, "*", "")
dimnames(sig_label) <- dimnames(mat_rho)

# ---------------------------------------------------------
# ---------------------------------------------------------
new_names <- mag2label[rownames(mat_rho)]
new_names[is.na(new_names)] <- rownames(mat_rho)
new_names <- make.unique(new_names)

rownames(mat_rho)   <- new_names
rownames(sig_label) <- new_names
rownames(row_anno)  <- new_names

# ---------------------------------------------------------
# ---------------------------------------------------------
bk   <- seq(-1, 1, length = 101)
cols <- colorRampPalette(c("blue", "white", "red"))(100)
row_anno$Traits <- recode(row_anno$Traits,
                            "Rib_eye" = "Ribeye area",
                            "Weight" = "Carcass weight")


tiff(
  "Heatmap_WeightMAGs_VFA_PosNeg_withGenus.tiff",
  width = 1800,
  height = 2500,
  res = 300,
  compression = "lzw"
)

pheatmap(
  mat_rho,
  color = cols,
  breaks = bk,
  display_numbers = sig_label,   # ← * 表示
  number_color = "black",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_anno,
  annotation_colors = ann_cols,
  fontsize_row = 8,
  fontsize_col = 10,
  main = ""
)

dev.off()
