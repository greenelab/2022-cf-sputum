## ---------------------------
##
## Script name: snakemake_combine_compendia.R
##
## Description: Combine TPM and NumReads from new samples with existing raw compendia
##
## Snakemake inputs: csv files for NumReads and TPM from new samples
##                   csv files for raw NumReads and TPM for existing Pa compendia
##
## Author: Taylor Reiter
##
## Date Created: 2022-03-11
##
## ---------------------------

library(dplyr)
library(readr)
library(janitor)
# read in raw original compendia ------------------------------------
tpm_og <- read_csv(snakemake@input[["tpm_og"]], show_col_types = FALSE) %>%
# tpm_og <- read_csv("inputs/original_compendia/TPM_pa14_cdna_k15.csv") %>%
  rename(TxName = "...1") 
colnames(tpm_og) <- gsub("\\.salmon", "", basename(colnames(tpm_og))) # fix colnames so only SRX
tpm_og <- tpm_og %>%
  clean_names(case = "all_caps")


numreads_og <- read_csv(snakemake@input[["numreads_og"]], show_col_types = FALSE) %>%
# numreads_og <- read_csv("inputs/original_compendia/num_reads_pao1_cdna_k15.csv") %>%
  rename(TxName = "...1")
colnames(numreads_og) <- gsub("\\.salmon", "", basename(colnames(numreads_og)))  # fix colnames so only SRX
numreads_og <- numreads_og %>%
  clean_names(case = "all_caps")

# read in raw new srx samples ---------------------------------------
tpm_new <- read_csv(snakemake@input[["tpm_new"]], show_col_types = FALSE) %>%
#tpm_new <- read_csv("outputs/combined_new_srx/TPM_pa14.csv") %>%
  clean_names(case = "all_caps")
numreads_new <- read_csv(snakemake@input[["numreads_new"]], show_col_types = FALSE) %>%
# numreads_new <- read_csv("outputs/combined_new_srx/num_reads_pao1.csv") %>%
  clean_names(case = "all_caps")

# check that no new samples are already in the compendia ------------------

new_samples_in_og_compendia <- colnames(numreads_new)[colnames(numreads_new[2:ncol(numreads_new)]) %in% colnames(numreads_og[3:ncol(numreads_og)])]
print(paste0(length(new_samples_in_og_compendia), " new samples were already in the original compendia."))

# some of the new samples are already in the compendia, remove them before joining
if(length(new_samples_in_og_compendia) > 0){
  numreads_new <- numreads_new %>%
    select(-one_of(colnames(numreads_og[3:ncol(numreads_og)])))
  
  tpm_new <- tpm_new %>%
    select(-one_of(colnames(tpm_og[3:ncol(tpm_og)])))
  
  print(paste0("Filtered ", length(new_samples_in_og_compendia), " new samples that were already in the original compendia."))
}
# combine and write out as csv --------------------------------------------

# use left join on old compendia gene names; will filter any additional names
# that sneak into the new compendia, although I can't think of a time this would
# actually happen
tpm <- left_join(tpm_og, tpm_new, by = "NAME") %>%
  select(TxName = TX_NAME, Name = NAME, starts_with("E"), starts_with("S"))
write_csv(tpm, snakemake@output[['tpm']])
numreads <- left_join(numreads_og, numreads_new, by = "NAME") %>%
  select(TxName = TX_NAME, Name = NAME, starts_with("E"), starts_with("S"))
write_csv(numreads, snakemake@output[['numreads']])