library(readr)
library(dplyr)
library(tibble)

# read in counts and metadata ---------------------------------------------

# set strain using wildcard so that proper set of control samples can be extracted
strain <- unlist(snakemake@wildcards[['strain']])

# read in hogan lab sample counts
#counts <- read_csv("outputs/combined_new/num_reads_pao1.csv")
counts <- read_csv(snakemake@input[['counts']], show_col_types = F)

# read in hogan lab sample metadata
#metadata <- read_csv("inputs/hogan_metadata.csv")
metadata <- read_csv(snakemake@input[['metadata']], show_col_types = F)

# determine experiment that is being compared ----------------------------------

# figure out what comparison should be formatted with this run by parsing the 
# snakemake wildcard. Each "hogan_comparison" contains two experiments, separated
# by "-vs-".
comparison <- unlist(snakemake@wildcards[['hogan_comparison']])
# comparison <- 'asm-vs-asm_m'

# parse comparison to determine which two sample sets to place in template files
set1 <- gsub("-.*", "", comparison)
set2 <- gsub(".*-vs-", "", comparison)

# separate hogan lab samples into sets based on sample type --------------------

spu_m_samples <- metadata %>%
  filter(metals == "metals") %>%
  filter(!grepl("ASM", x = sample))

asm_m_samples <- metadata %>%
  filter(metals == "metals") %>%
  filter(grepl("ASM", x = sample))

spu_samples <- metadata %>%
  filter(metals == "no metals") %>%
  filter(!grepl("ASM", x = sample)) %>%
  filter(!grepl("^M", x = sample))

asm_samples <- metadata %>%
  filter(metals == "no metals") %>%
  filter(grepl("ASM", x = sample))

m63_samples <- metadata %>%
  filter(metals == "no metals") %>%
  filter(grepl("^M", x = sample))

# spu_m_counts <- counts %>%
#   select(Name, one_of(metal_samples$sample)) %>% # select Name and metal samples
#   left_join(compendium, by = "Name") %>% # join to compendium controls
#   column_to_rownames("Name") %>%       # put gene name as rowname, so it will become colname
#   t()                                  # transpose

spu_m_counts <- counts %>% select(Name, one_of(spu_m_samples$sample)) 

asm_m_counts <- counts %>% select(Name, one_of(asm_m_samples$sample)) # select Name and metal samples

spu_counts <- counts %>% select(Name, one_of(spu_samples$sample)) # select Name and nonmetal samples
  
asm_counts <- counts %>% select(Name, one_of(asm_samples$sample)) # select Name and artificial samples

m63_counts <- counts %>% select(Name, one_of(m63_samples$sample)) # select Name and m63 samples


# create metadata and count tables with correct contrasts ----------------------

if(comparison == 'asm-vs-asm_m'){
  counts_out <- left_join(asm_counts, asm_m_counts, by = "Name") %>% # join by gene name
    column_to_rownames("Name") %>%       # put gene name as rowname, so it will become colname
    t()                                  # transpose
  groups_out <- data.frame(sample = rownames(counts_out)) %>%
    mutate(group = ifelse(sample %in% colnames(asm_counts), 1, 2))
} else if(comparison == "spu-vs-spu_m"){
  counts_out <- left_join(spu_counts, spu_m_counts, by = "Name") %>% # join by gene name
    column_to_rownames("Name") %>%       # put gene name as rowname, so it will become colname
    t()                                  # transpose
  groups_out <- data.frame(sample = rownames(counts_out)) %>%
    mutate(group = ifelse(sample %in% colnames(spu_counts), 1, 2))
} else if(comparison == "spu-vs-asm"){
  counts_out <- left_join(spu_counts, asm_counts, by = "Name") %>% # join by gene name
    column_to_rownames("Name") %>%       # put gene name as rowname, so it will become colname
    t()                                  # transpose
  groups_out <- data.frame(sample = rownames(counts_out)) %>%
    mutate(group = ifelse(sample %in% colnames(spu_counts), 1, 2))
} else if(comparison == "spu-vs-m63"){
  counts_out <- left_join(spu_counts, m63_counts, by = "Name") %>% # join by gene name
    column_to_rownames("Name") %>%       # put gene name as rowname, so it will become colname
    t()                                  # transpose
  groups_out <- data.frame(sample = rownames(counts_out)) %>%
    mutate(group = ifelse(sample %in% colnames(spu_counts), 1, 2))
} else if(comparison == "spu_m-vs-asm_m"){
  counts_out <- left_join(spu_m_counts, asm_m_counts, by = "Name") %>% # join by gene name
    column_to_rownames("Name") %>%       # put gene name as rowname, so it will become colname
    t()                                  # transpose
  groups_out <- data.frame(sample = rownames(counts_out)) %>%
    mutate(group = ifelse(sample %in% colnames(spu_m_counts), 1, 2))

} else if(comparison == "asm-vs-m63"){
  counts_out <- left_join(asm_counts, m63_counts, by = "Name") %>% # join by gene name
    column_to_rownames("Name") %>%       # put gene name as rowname, so it will become colname
    t()                                  # transpose
  groups_out <- data.frame(sample = rownames(counts_out)) %>%
    mutate(group = ifelse(sample %in% colnames(asm_counts), 1, 2))
}

# ponyo_out depends on groups_out and the comparison, so it can be run outside
# of the if/else statements.
ponyo_out  <- data.frame(experiment_colname = paste0(strain, "_", comparison), 
                         sample_id_colname  = groups_out$sample)
# write results -----------------------------------------------------------

# use base R write.table() to write counts so that rownames are written as an index
write.table(counts_out, snakemake@output[['num-reads']], quote = F, sep = "\t")

write_tsv(groups_out, snakemake@output[['groups']])

write_csv(ponyo_out, snakemake@output[['ponyo']])
