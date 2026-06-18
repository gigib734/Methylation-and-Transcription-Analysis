methylation_data <- read.csv("C:/Users/gbloc/OneDrive/Desktop/PCDHG_Gila/Coding/Methylation & Transcription Correlation/Meta-analysis results DMPs blood PCDH genes - Meta-analysis results DMPs blood PCDH genes (1).csv")
transcription_data <- read.csv("C:/Users/gbloc/OneDrive/Desktop/PCDHG_Gila/Coding/Methylation & Transcription Correlation/Meta-analysis results DEGs blood PCDH genes - Meta-analysis results DEGs blood PCDH genes.csv.csv")
Filtered_methylation_data <- methylation_data %>%
  select(CpG, Annotated_genes, chr_state, Effect_size, SE, CGI_position)
library(dplyr)
library(tidyverse)
library(readr)
library(ggplot2)

#separating annotated genes out
Gene_Filtered_Methylation_data <- Filtered_methylation_data %>%
  separate_rows(Annotated_genes, sep = ";") %>%
  mutate(Annotated_genes = str_trim(Annotated_genes))

#Filtering by cpg region and chromatin state


#meaneffectbygene
Mean_gene_effect <- Gene_Filtered_Methylation_data %>%
  group_by(Annotated_genes) %>%
  summarise(
    methylation_mean = mean(Effect_size, na.rm = TRUE),
    methylation_median = median(Effect_size, na.rm = TRUE),
    n_CpG = n(),
    Mean_Effect_Size = mean(Effect_size, na.rm = TRUE),
  )

Gene_Filtered_Methylation_data <- Gene_Filtered_Methylation_data %>%
  mutate(
    Effect_size = as.numeric(Effect_size),
    SE = as.numeric(SE))

transcription_data<- transcription_data %>%
  rename(Annotated_genes = MarkerName)

dat <- inner_join(Mean_gene_effect, transcription_data, by = "Annotated_genes")
cor.test(dat$methylation_mean, dat$Effect, method = "spearman")

plot_data <- Mean_gene_effect %>%
  rename(gene = Annotated_genes) %>%
  inner_join(
    transcription_data %>% rename(gene = Annotated_genes),
    by = "gene"
  ) %>%
  select(gene, methylation_mean, Effect)
  )

head(plot_data)
ggplot(plot_data, aes(x = methylation_mean, y = Effect)) +
  geom_point(size = 3, color = "steelblue", alpha = 0.6) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(
    title = "Methylation Mean vs Transcription Effect",
    x = "Methylation Mean",
    y = "Transcription Effect"
  )

#weighted mean based on SE
Weighted_gene_effect <- Gene_Filtered_Methylation_data %>%
  filter(!is.na(SE) & SE > 0) %>%  # Remove invalid SE values
  group_by(Annotated_genes) %>%
  summarise(
    n_CpG = n(),
    methylation_mean = mean(Effect_size, na.rm = TRUE),
    weighted_mean = sum(Effect_size / SE^2, na.rm = TRUE) / sum(1 / SE^2, na.rm = TRUE),
    se_weighted = 1 / sqrt(sum(1 / SE^2, na.rm = TRUE)),
    .groups = 'drop'
  ) %>%
  arrange(desc(abs(weighted_mean)))

Hm <- inner_join(Weighted_gene_effect, transcription_data, by = "Annotated_genes")
cor.test(Hm$weighted_mean, dat$Effect, method = "spearman")

plot_weighted_data <- Weighted_gene_effect %>%
  rename(gene = Annotated_genes) %>%
  inner_join(
    transcription_data %>% rename(gene = Annotated_genes),
    by = "gene"
  ) %>%
  select(gene, weighted_mean, Effect)


ggplot(plot_weighted_data, aes(x = weighted_mean, y = Effect)) +
  geom_point(size = 3, color = "steelblue", alpha = 0.6) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(
    title = "Weighted_Methylation Mean vs Transcription Effect",
    x = "Methylation Mean",
    y = "Transcription Effect"
  )
#Filtering of Methylation data only using TSS regions and island and island shores

Island_Filtered_methylation_data <- methylation_data %>%
  filter(str_detect(chr_state,"(?i)TSS"))
head(Region_filtered_data)
nrow(Island_Filtered_methylation_data)

Island_sun$Effect_size <- as.numeric(Island_sun$Effect_size)
str(IIsland_Filtered_methylation_datastr(Island_sun$Effect_size)

Island_sun <- Island_Filtered_methylation_data %>%
  filter(CGI_position %in% c("Island", "N_Shore", "S_Shore"))
nrow(Island_sun)

Island_sun <- Island_sun %>%
  mutate(Effect_size = as.numeric(Effect_size))

Island_sun <- Island_sun %>%
  separate_rows(Annotated_genes, sep = ";") %>%
  mutate(Annotated_genes = str_trim(Annotated_genes))


Island_mean <- Island_sun %>%
  group_by(Annotated_genes) %>%
  summarise(
    methylation_mean = mean(Effect_size, na.rm = TRUE),
    n_CpG = n(),
    .groups = 'drop'
  )

Island_weighted_mean <- Island_sun %>%
  filter(!is.na(SE) & SE > 0) %>%  # Remove invalid SE values
  group_by(Annotated_genes) %>%
  summarise(
    n_CpG = n(),
    methylation_mean = mean(Effect_size, na.rm = TRUE),
    weighted_mean = sum(Effect_size / SE^2, na.rm = TRUE) / sum(1 / SE^2, na.rm = TRUE),
    se_weighted = 1 / sqrt(sum(1 / SE^2, na.rm = TRUE)),
    .groups = 'drop'
  ) %>%
  arrange(desc(abs(weighted_mean)))


island_1 <- inner_join(Island_mean, transcription_data, by = "Annotated_genes")
cor.test(island_1$methylation_mean, island_1$Effect, method = "spearman")

#plot the unweighted island meth vs trans
ggplot(island_1, aes(x = methylation_mean, y = Effect)) +
  geom_point(size = 3, color = "steelblue", alpha = 0.6) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(
    title = "Island TSS Methylation Mean vs Transcription Effect",
    x = "Methylation Mean",
    y = "Transcription Effect"
  )

ggplot(island_1, aes(x = methylation_mean, y = Effect)) +
  geom_point(size = 3, color = "steelblue", alpha = 0.6) +
  geom_smooth(method = "lm") +
  geom_text(
    aes(label = ifelse(
      Annotated_genes %in% c(
        "PCDHGA1",
        "PCDHGA2",
        "PCDHGA3",
        "PCDHGA4",
        "PCDHGB1"
      ),
      Annotated_genes,
      ""
    )),
    vjust = -0.5
  ) +
  theme_minimal() +
  labs(
    title = "Island TSS Methylation Mean vs Transcription Effect",
    x = "Methylation Mean",
    y = "Transcription Effect"
  )


#plot the weightedmmeth vs trans fr island and ts filter
island_2 <- inner_join(Island_weighted_mean, transcription_data, by = "Annotated_genes")
cor.test(island_2$weighted_mean, island_2$Effect, method = "spearman")


ggplot(island_2, aes(x = methylation_mean, y = Effect)) +
  geom_point(size = 3, color = "steelblue", alpha = 0.6) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(
    title = "Island TSS Weighted Methylation Mean vs Transcription Effect",
    x = "Methylation Mean",
    y = "Transcription Effect"
  )

unique(methylation_data$CGI_position)
unique(methylation_data$chr_state)

#mean and weighted mean when just filtering for pcdg target genes

PCDH_target <- Gene_Filtered_Methylation_data %>%
   filter(Annotated_genes %in% c("PCDHGA1","PCDHGA2","PCDHGA3","PCDHGA4","PCDHGB1"))

unique(PCDH_target$Annotated_genes)

nrow(PCDH_target)


PCDH_by_gene <- PCDH_target %>%
  group_by(Annotated_genes) %>%
  summarise(
    methylation_mean = mean(Effect_size, na.rm = TRUE),
    methylation_median = median(Effect_size, na.rm = TRUE),
    n_CpG = n(),
    Mean_Effect_Size = mean(Effect_size, na.rm = TRUE),
  )

PCDH_mean <- inner_join(PCDH_by_gene, transcription_data, by = "Annotated_genes")
cor.test(PCDH_mean$methylation_mean, PCDH_mean$Effect, method = "spearman")

ggplot(PCDH_mean, aes(x = methylation_mean, y = Effect)) +
  geom_point(size = 3, color = "steelblue", alpha = 0.6) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(
    title = "PCDHG Target Methylation Mean vs Transcription Effect",
    x = "Methylation Mean",
    y = "Transcription Effect"
  )

PCDH_weighted <- PCDH_target %>%
  filter(!is.na(SE) & SE > 0) %>%  # Remove invalid SE values
  group_by(Annotated_genes) %>%
  summarise(
    n_CpG = n(),
    methylation_mean = mean(Effect_size, na.rm = TRUE),
    weighted_mean = sum(Effect_size / SE^2, na.rm = TRUE) / sum(1 / SE^2, na.rm = TRUE),
    se_weighted = 1 / sqrt(sum(1 / SE^2, na.rm = TRUE)),
    .groups = 'drop'
  ) %>%
  arrange(desc(abs(weighted_mean)))

PCDH_weighted_mean <- inner_join(PCDH_weighted, transcription_data, by = "Annotated_genes")
cor.test(PCDH_weighted_mean$weighted_mean, PCDH_weighted_mean$Effect, method = "spearman")

ggplot(PCDH_weighted_mean, aes(x = weighted_mean, y = Effect)) +
  geom_point(size = 3, color = "steelblue", alpha = 0.6) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(
    title = "PCDHG Target Weighted Methylation Mean vs Transcription Effect",
    x = "Methylation Mean",
    y = "Transcription Effect"
  )


#can i take my target genes an do the same filtering for TSS and island regions?

Island_sun_pcdh <- Island_sun %>%
  filter(
    Annotated_genes %in% c(
      "PCDHGA1",
      "PCDHGA2",
      "PCDHGA3",
      "PCDHGA4",
      "PCDHGB1"
    )
  )


Island_sun_pcdh <- Island_sun_pcdh %>%
  separate_rows(Annotated_genes, sep = ";") %>%
  mutate(Annotated_genes = str_trim(Annotated_genes))


Island_pcdh_mean <- Island_sun_pcdh %>%
  group_by(Annotated_genes) %>%
  summarise(
    methylation_mean = mean(Effect_size, na.rm = TRUE),
    n_CpG = n(),
    .groups = 'drop'
  )

PCDH_island_mean_thing <- inner_join(Island_pcdh_mean, transcription_data, by = "Annotated_genes")
cor.test(PCDH_island_mean_thing$methylation_mean, PCDH_island_mean_thing$Effect, method = "spearman")

ggplot(PCDH_island_mean_thing, aes(x = methylation_mean, y = Effect)) +
  geom_point(size = 3, color = "steelblue", alpha = 0.6) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(
    title = "PCDHG Target Methylation Mean vs Transcription Effect",
    x = "Methylation Mean",
    y = "Transcription Effect"
  )
