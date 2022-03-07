library(dplyr)
library(readr)
library(stringr)
setwd("~/github/2022-cf-sputum/")

pa14 <- read_csv("inputs/pa_compendia/Dataset_S4_PA14_compendium.gz")

accession_in_comp_all <- data.frame()
for(accession in metadata$experiment_accession){
  accession_in_comp <- str_detect(string = colnames(pa14), pattern = accession)
  accession_in_comp <- TRUE %in% accession_in_comp
  accession_in_comp_df <- data.frame(accession, accession_in_comp)
  accession_in_comp_all <- bind_rows(accession_in_comp_all, accession_in_comp_df)
}

srx <- accession_in_comp_all %>%
  filter(accession_in_comp == F)
cat(srx$accession, sep = "', '")
