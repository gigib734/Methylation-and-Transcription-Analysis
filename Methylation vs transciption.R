methylation_data <- read.csv("C:/Users/gbloc/OneDrive/Desktop/PCDHG_Gila/Coding/Methylation & Transcription Correlation/Meta-analysis results DMPs blood PCDH genes - Meta-analysis results DMPs blood PCDH genes (1).csv")
transcription_data <- read.csv("C:/Users/gbloc/OneDrive/Desktop/PCDHG_Gila/Coding/Methylation & Transcription Correlation/Meta-analysis results DEGs blood PCDH genes - Meta-analysis results DEGs blood PCDH genes.csv.csv")
Filtered_methylation_data <- methylation_data %>%
  select(CpG, Annotated_genes, chr_state, Effect_size)
library(dplyr)
library(tidyverse)
library(readr)

#Filtering of Methylation data to do it by gene

str(Filtered_methylation_data)
#separating annotated genes out
Gene_Filtered_Methylation_data <- Filtered_methylation_data %>%
  separate_rows(Annotated_genes, sep = ";") %>%
  mutate(Annotated_genes = str_trim(Annotated_genes))


#Filtering of Methylation data only using TSS regions


