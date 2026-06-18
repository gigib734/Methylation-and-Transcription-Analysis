methylation_data <- read.csv("C:/Users/gbloc/OneDrive/Desktop/PCDHG_Gila/Coding/Methylation & Transcription Correlation/Meta-analysis results DMPs blood PCDH genes - Meta-analysis results DMPs blood PCDH genes (1).csv")
transcription_data <- read.csv("C:/Users/gbloc/OneDrive/Desktop/PCDHG_Gila/Coding/Methylation & Transcription Correlation/Meta-analysis results DEGs blood PCDH genes - Meta-analysis results DEGs blood PCDH genes.csv.csv")
Filtered_methylation_data <- methylation_data %>%
  select(CpG, Annotated_genes, chr_state, Effect_size, SE)
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

Region_filtered_data <- methylation_data %>%
  filter(str_detect(chr_state,"(?i)TSS"))
head(Region_filtered_data)
nrow(Region_filtered_data)

Island_sun <- Region_filtered_data %>%
  filter(CGI_position %in% c("Island", "N_Shore", "S_Shore"))
nrow(Island_sun)



unique(methylation_data$CGI_position)
unique(methylation_data$chr_state)
