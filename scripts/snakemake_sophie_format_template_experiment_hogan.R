library(readr)
library(dplyr)
library(tibble)

# read in counts and metadata ---------------------------------------------

# set strain using wildcard so that proper set of control samples can be extracted
strain <- unlist(snakemake@wildcards[['strain']])

# read in hogan lab sample counts
#counts <- read_csv("outputs/combined_new/num_reads_pa14.csv")
counts <- read_csv(snakemake@input[['counts']], show_col_types = F)

# read in hogan lab sample metadata
#metadata <- read_csv("inputs/hogan_metadata.csv")
metadata <- read_csv(snakemake@input[['metadata']], show_col_types = F)

# read in compendium counts so that "control" counts can be extracted
#compendium <- read_csv(snakemake@input[['compendium']], show_col_types = F) %>%
compendium <- read_csv("inputs/original_compendia/num_reads_pao1_cdna_k15.csv", show_col_types = F) %>%
  rename(tx_name = "...1") %>% # fix auto named col
  select(-tx_name) %>%  # remove txname, as its not needed
  rename_with(~ gsub("\\.salmon", "", basename(.x))) # edit colnames to exp accession


# determine control set of samples and extract from compendium ------------

# control samples were selected based on the following criteria:
# 1. same strain as compendium (i.e. either PAO1 or PA14)
# 2. WT cells grown in LB medium
# 3. Sequenced in biological triplicate

if(strain == "pao1") {
  # If PAO1, use WT from PRJNA613826 as controls for differential expression
  controls <- c("SRX8487076", "SRX8487077", "SRX8487078")
  # grab just the controls out of the compendium
  compendium <- compendium %>%
    select(Name, one_of(controls))
} else {
  # If PA14, use WT from PRJNA291596
  # alternative experiments I considered were:
  # 1. PRNA343895 (unclear if LB medium was used)
  # 2. SRX4641857, SRX4641856 (sequenced in duplicate)
  # 3. SRX8030428, SRX8030427, SRX8030426, SRX8030425 (SRX8030427 and SRX8030425 are tech reps)
  controls <- c("SRX1127445", "SRX1127446", "SRX1127447")
  # grab just the controls out of the compendium
  compendium <- compendium %>%
    select(Name, one_of(controls))
}

# separate hogan lab samples, metal nonmetals -----------------------------

metal_samples <- metadata %>%
  filter(metals == "metals") %>%
  filter(!grepl("ASM", x = sample))

artificial_metal_samples <- metadata %>%
  filter(metals == "metals") %>%
  filter(grepl("ASM", x = sample))

nonmetal_samples <- metadata %>%
  filter(metals == "no metals") %>%
  filter(!grepl("ASM", x = sample)) %>%
  filter(!grepl("^M", x = sample))

artificial_nonmetal_samples <- metadata %>%
  filter(metals == "no metals") %>%
  filter(grepl("ASM", x = sample))

m_nonmetal_samples <- metadata %>%
  filter(metals == "no metals") %>%
  filter(grepl("^M", x = sample))

counts_metals <- counts %>%
  select(Name, one_of(metal_samples$sample)) %>% # select Name and metal samples
  left_join(compendium, by = "Name") %>% # join to compendium controls
  column_to_rownames("Name") %>%       # put gene name as rowname, so it will become colname
  t()                                  # transpose

counts_artificial_metals <- counts %>%
  select(Name, one_of(artificial_metal_samples$sample)) %>% # select Name and metal samples
  left_join(compendium, by = "Name") %>% # join to compendium controls
  column_to_rownames("Name") %>%       # put gene name as rowname, so it will become colname
  t()                                  # transpose

counts_nonmetals <- counts %>%
  select(Name, one_of(nonmetal_samples$sample)) %>% # select Name and nonmetal samples
  left_join(compendium, by = "Name") %>% # join to compendium controls
  column_to_rownames("Name") %>%       # put gene name as rowname, so it will become colname
  t()                                  # transpose

counts_artificial_nonmetals <- counts %>%
  select(Name, one_of(artificial_nonmetal_samples$sample)) %>% # select Name and nonmetal samples
  left_join(compendium, by = "Name") %>% # join to compendium controls
  column_to_rownames("Name") %>%       # put gene name as rowname, so it will become colname
  t()                                  # transpose

counts_m_nonmetals <- counts %>%
  select(Name, one_of(m_nonmetal_samples$sample)) %>% # select Name and nonmetal samples
  left_join(compendium, by = "Name") %>% # join to compendium controls
  column_to_rownames("Name") %>%       # put gene name as rowname, so it will become colname
  t()                                  # transpose

# create metadata tables --------------------------------------------------

# create a data frame to record the sample groupings, which will be used by 
# sophie to perform differential expression analysis

sophie_metadata_metals <- data.frame(sample = rownames(counts_metals)) %>%
  mutate(group = ifelse(sample %in% controls, 1, 2))

sophie_metadata_artificial_metals <- data.frame(sample = rownames(counts_artificial_metals)) %>%
  mutate(group = ifelse(sample %in% controls, 1, 2))

sophie_metadata_nonmetals <- data.frame(sample = rownames(counts_nonmetals)) %>%
  mutate(group = ifelse(sample %in% controls, 1, 2))

sophie_metadata_artificial_nonmetals <- data.frame(sample = rownames(counts_artificial_nonmetals)) %>%
  mutate(group = ifelse(sample %in% controls, 1, 2))

sophie_metadata_m_nonmetals <- data.frame(sample = rownames(counts_m_nonmetals)) %>%
  mutate(group = ifelse(sample %in% controls, 1, 2))


# create ponyo metadata tables --------------------------------------------

# experiment_colname | sample_id_colname
# pao1_hogan_metals | sample_id_0
# pao1_hogan_metals | sample_id_1
# pao1_hogan_metals | sample_id_2
# pao1_hogan_metals | sample_id_3

ponyo_metadata_metals <- data.frame(experiment_colname = "pao1_hogan_metals", 
                                    sample_id_colname  = sophie_metadata_metals$sample)
write_csv(ponyo_metadata_metals, snakemake@output[['ponyo_metals']])
# write results -----------------------------------------------------------

# use base R write.table() to write counts so that rownames are written as an index
write.table(counts_metals, snakemake@output[['metals']], quote = F, sep = "\t")
write.table(counts_artificial_metals, snakemake@output[['artificial_metals']], quote = F, sep = "\t")
write.table(counts_nonmetals, snakemake@output[['sputum']], quote = F, sep = "\t")
write.table(counts_artificial_nonmetals, snakemake@output[['artificial_sputum']], quote = F, sep = "\t")
write.table(counts_m_nonmetals, snakemake@output[['m_sputum']], quote = F, sep = "\t")

write_tsv(sophie_metadata_metals, snakemake@output[['meta_metals']])
write_tsv(sophie_metadata_artificial_metals, snakemake@output[['meta_artificial_metals']])
write_tsv(sophie_metadata_nonmetals, snakemake@output[['meta_sputum']])
write_tsv(sophie_metadata_artificial_nonmetals, snakemake@output[['meta_artificial_sputum']])
write_tsv(sophie_metadata_m_nonmetals, snakemake@output[['meta_m_sputum']])
