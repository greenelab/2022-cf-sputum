library(readr)
library(dplyr)
library(tibble)

# read in counts and metadata ---------------------------------------------

# set strain using wildcard so that proper set of control samples can be extracted
strain <- unlist(snakemake@wildcards[['strain']])
# strain <- "pao1"

# read in hogan lab sample counts
#counts <- read_csv("outputs/combined_compendia/num_reads_pao1.csv")
counts <- read_csv(snakemake@input[['counts']], show_col_types = F)

# read in hogan lab sample metadata
#metadata <- read_csv("inputs/metadata.csv")
metadata <- read_csv(snakemake@input[['metadata']], show_col_types = F) %>%
  filter(!experiment_accession %in% c("SRX7101177", "SRX7101178", "SRX7101179", "SRX7101180", "SRX7101181")) # filter out low count samples; see comment below

# below I recorded the counts for these samples when mapped against PAO1.
# there aren't enough counts so I removed them here.
# SRX7101177 SRX7101178 SRX7101179 SRX7101180 SRX7101181 
# 131        476         49        164        263
# Every other library has at least 7000 reads mapped to pao1.

# determine experiment that is being compared ----------------------------------

# figure out what comparison should be formatted with this run by parsing the 
# snakemake wildcard. Each "hogan_comparison" contains two experiments, separated
# by "-vs-".
comparison <- unlist(snakemake@wildcards[['pub_comparison']])
# comparison <- 'pub-vs-lbpao1'

# parse comparison to determine which two sample sets to place in template files
set1 <- gsub("-.*", "", comparison)
set2 <- gsub(".*-vs-", "", comparison)

# separate hogan lab samples into sets based on sample type --------------------

pub_samples <- metadata %>%
  filter(pa_in_reads == TRUE) 

m63_samples <- c("M631", "M632", "M633")

lbpao1_samples <- c('SRX5661364', 'SRX5661363', 'SRX6976948', 'SRX6976949',
                    'SRX6976950', 'SRX8862509', 'SRX8862510', 'ERX2813653',
                    'ERX2813654', 'ERX2813655', 'ERX2068559', 'ERX2068560',
                    'ERX2068561', 'SRX4579961', 'SRX4579962', 'SRX4579963',
                    'SRX8487076', 'SRX8487077', 'SRX8487078')

lbpa14_samples <-  c('SRX474156', 'SRX474157', 'SRX8030422', 'SRX8030424',
                     'SRX474130', 'SRX474131', 'SRX7299397', 'SRX7299398',
                     'SRX8030426', 'SRX8030428', 'SRX470383', 'SRX470384',
                     'SRX804118', 'SRX804119', 'SRX4641857', 'SRX4641858')


pub_counts <- counts %>% select(Name, one_of(pub_samples$experiment_accession)) 

lbpao1_counts <- counts %>% select(Name, one_of(lbpao1_samples)) 
  
lbpa14_counts <- counts %>% select(Name, one_of(lbpa14_samples)) 

m63_counts <- counts %>% select(Name, one_of(m63_samples)) # select Name and m63 samples


# create metadata and count tables with correct contrasts ----------------------

create_count_and_group_metadata <- function(counts1, counts2){
  # function that takes in two counts dfs, joins them into a single counts df, 
  # and uses sample names from the first counts dataframe to generate a groups dataframe. 
  # The groups dataframe is used by sophie/DESeq2 to perform differential expression analysis.
  # The samples in the first dataframe will be labelled as "1", 
  # and the samples in the second dataframe will be labelled as "2".
  counts_out <- left_join(counts1, counts2, by = "Name") %>% # join by gene name
    column_to_rownames("Name") %>%       # put gene name as rowname, so it will become colname
    t()                                  # transpose
  groups_out <- data.frame(sample = rownames(counts_out)) %>%
    mutate(group = ifelse(sample %in% colnames(counts1), 1, 2))
  return(list(counts = counts_out, groups = groups_out))
}

if(comparison == 'pub-vs-m63'){
  out <- create_count_and_group_metadata(pub_counts, m63_counts)
} else if(comparison == "pub-vs-lbpao1"){
  out <- create_count_and_group_metadata(pub_counts, lbpao1_counts)
} else if(comparison == "pub-vs-lbpa14"){
  out <- create_count_and_group_metadata(pub_counts, lbpa14_counts)
} 

# ponyo_out depends on groups_out and the comparison, so it can be run outside
# of the if/else statements.
ponyo_out  <- data.frame(experiment_colname = paste0(strain, "_", comparison), 
                         sample_id_colname  = out$groups$sample)
# write results -----------------------------------------------------------

# use base R write.table() to write counts so that rownames are written as an index
write.table(out$counts, snakemake@output[['num_reads']], quote = F, sep = "\t")

write_tsv(out$groups, snakemake@output[['grps']])

write_csv(ponyo_out, snakemake@output[['ponyo']])
