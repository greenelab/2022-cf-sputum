## ---------------------------
##
## Script name: snakemake_quant_collect.R
##
## Description: Combines all quant.sf files into a csv as part of a snakemake pipeline
##
## Snakemake inputs: all quant.sf files for one strain (strain controlled by snakemake)
##                   annotation map from transcript ID to gene name
##                   output csv file names for NumReads and TPM
##
## Author: Taylor Reiter
##
## Date Created: 2022-03-11
##
## ---------------------------

library(dplyr)
library(readr)
library(purrr)
library(tidyr)

annotation_map <- read_csv(snakemake@input[["annot_map"]], col_names = c("tx_name", "tx_name2", "Name")) %>%
  select(-tx_name2)
#annotation_map <- read_csv("~/Downloads/pa14_gene_names.csv", col_names = c("tx_name", "tx_name2", "Name")) %>%
#  select(-tx_name2)
sf_files <- unlist(snakemake@input[["quant"]])
#sf_files <- Sys.glob("outputs/salmon/pa14/*/quant.sf")

sf <- sf_files %>%
  set_names() %>%                                                    # add filename as a column
  map_dfr(read_tsv, show_col_types = FALSE, .id = "srx") %>%                                 # read in files
  mutate(srx = basename(dirname(srx))) %>%                           # edit filename to srx ID
  left_join(annotation_map, by = c("Name" = "tx_name")) %>%          # join to gene name
  select(srx, Name = Name.y, Length, EffectiveLength, TPM, NumReads) # drop transcript name

numreads <- sf %>%
  select(srx, Name, NumReads) %>%                                       # select NumReads
  pivot_wider(id_cols = Name, values_from = NumReads, names_from = srx) # transform long to wide
  
tpm <- sf %>%
  select(srx, Name, TPM) %>%                                       # select TPM
  pivot_wider(id_cols = Name, values_from = TPM, names_from = srx) # transform long to wide

# write out tables for NumReads and TPM
write_csv(numreads, snakemake@output[['numreads']])
write_csv(tpm, snakemake@output[['tpm']])