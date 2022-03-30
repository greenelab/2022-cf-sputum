library(readr)
library(dplyr)


# read in counts and metadata ---------------------------------------------

#counts <- read_csv("outputs/combined_new/num_reads_pa14.csv")
counts <- read_csv(snakemake@input[['counts']], show_col_types = F)

#metadata <- read_csv("inputs/hogan_metadata.csv")
metadata <- read_csv(snakemake@input[['metadata']], show_col_types = F)


# separate hogan lab samples, metal nonmetals -----------------------------

metal_samples <- metadata %>%
  filter(metals == "metals")

nonmetal_samples <- metadata %>%
  filter(metals == "no metals")

counts_metals <- counts %>%
  select(Name, one_of(metal_samples$sample)) %>%
  column_to_rownames("Name") %>%
  t()

counts_nonmetals <- counts %>%
  select(Name, one_of(nonmetal_samples$sample)) %>%
  column_to_rownames("Name") %>%
  t()


# write results -----------------------------------------------------------

# use base R write.csv() so that rownames are written as an index
write.csv(counts_metals, snakemake@output[['metals']], quote = F)
write.csv(counts_nonmetals, snakemake@output[['sputum']], quote = F)