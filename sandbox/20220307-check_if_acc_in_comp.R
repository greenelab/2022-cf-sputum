library(dplyr)
library(readr)
library(stringr)
setwd("~/github/2022-cf-sputum/")

metadata <- read_csv("inputs/metadata.csv")
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

metadata <- metadata %>%
  left_join(accession_in_comp_all, by = c("experiment_accession" = "accession"))

write_csv(metadata, "inputs/metadata.csv")

cat(srx$accession, sep = "', '")
