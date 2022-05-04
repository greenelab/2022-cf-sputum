library(readr)
library(dplyr)
library(tibble)

# this script reads in a Pa compendium and SRA-derived strain annotations for
# samples in the compendium. It filters the compendium to strain-specific samples
# and formats it for input to SOPHIE.

# read in compendium and metadata -----------------------------------------

strain <- unlist(snakemake@wildcards[['strain']])
strain <- toupper(strain) # fix case to match whats recorded in strain metadata tbl
#strain <- "pa14"

strain_metadata <- read_tsv(snakemake@input[['strain']], show_col_types = F) %>%
# strain_metadata <- read_tsv("inputs/original_compendia/SRA_annotations.tsv", show_col_types = F) %>%
  rename(experiment = Experiment, strain_type = "Strain type")

compendium <- read_csv(snakemake@input[['compendium']], show_col_types = F) %>%
# compendium <- read_csv("inputs/original_compendia/num_reads_pa14_cdna_k15.csv", show_col_types = F) %>%
  rename(tx_name = "...1") %>% # fix auto named col
  select(-tx_name) %>%  # remove txname, as its not needed
  rename_with(~ gsub("\\.salmon", "", basename(.x))) # edit colnames to exp accession

# subset to strain of interest ---------------------------------------------

# filter to curated strain type
strain_metadata_subset <- strain_metadata %>%
  filter(strain_type %in% strain)

# subset compendium to curated set; format for sophie (genes as colnames, samples as rownames)
compendium_subset <- compendium %>%
  select(Name, one_of(strain_metadata_subset$experiment)) %>%
  column_to_rownames("Name") %>%
  t()

# use base R write.csv() so rownames will be written w/o a colname, which I think
# will play best with sophie/python
write.table(compendium_subset, snakemake@output[['compendium']], quote = F, sep = "\t")
