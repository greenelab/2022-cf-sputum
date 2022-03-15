## ---------------------------
##
## Script name: snakemake_filter_and_normalize_compendia.R
##
## Description: Filter and normalize the compendia
## This script was modified from https://github.com/georgiadoing/pa-seq-compendia/blob/main/qc_filtering/pa-seq-compendia-QC-apply.Rmd
## to work with snakemake and tidyverse
##
## Snakemake inputs: csv files for NumReads and TPM that represent the combined new and original compendia
##                   filter functions
##
## Author: Taylor Reiter
##
## Date Created: 2022-03-11
##
## ---------------------------
library(readr)
library(dplyr)
library(tidyr)
source("scripts/filter_functions.R")
source("scripts/snakemake_annotation_functions.R")

# read in raw combined compendia files ------------------------------------

rnaseq_pao1_tpm <- read.csv(snakemake@input[['tpm']][[1]],
#rnaseq_pao1_tpm <- read.csv("outputs/combined_compendia/TPM_pao1.csv", 
                            stringsAsFactors = F, row.names = 1)[,-1]
rownames(rnaseq_pao1_tpm) <- make.unique(sapply(rownames(rnaseq_pao1_tpm), 
                                                function(x) cDNA_to_PAO1(x)))
rnaseq_pao1_counts <- read.csv(snakemake@input[['numreads']][[1]], 
#rnaseq_pao1_counts <- read.csv('outputs/combined_compendia/num_reads_pao1.csv', 
                               stringsAsFactors = F, row.names = 1)[,-1]
rownames(rnaseq_pao1_counts) <- make.unique(sapply(rownames(rnaseq_pao1_counts), 
                                                   function(x) cDNA_to_PAO1(x)))
# PA14
rnaseq_pa14_tpm <- read.csv(snakemake@input[['tpm']][[2]], 
#rnaseq_pa14_tpm <- read.csv('outputs/combined_compendia/tpm_pa14.csv', 
                            stringsAsFactors = F, row.names = 1)[,-1]
rownames(rnaseq_pa14_tpm) <- make.unique(sapply(rownames(rnaseq_pa14_tpm), 
                                                function(x) cDNA_to_PA14(x)))
rnaseq_pa14_counts <- read.csv(snakemake@input[['tpm']][[2]],
#rnaseq_pa14_counts <- read.csv('outputs/combined_compendia/num_reads_pa14.csv', 
                               stringsAsFactors = F, row.names = 1)[,-1]
rownames(rnaseq_pa14_counts) <- make.unique(sapply(rownames(rnaseq_pa14_counts), 
                                                   function(x) cDNA_to_PA14(x)))


# log ---------------------------------------------------------------------

# For filtering applications, use log-transformed data since it resembles more
# normal distributions. the hlog ("happy log") function takes the log of a 
# dataframe making sure the resulting data frame doesn't have any -Inf or NaN values.

rnaseq_pao1_tpm_log <- hlog(rnaseq_pao1_tpm)
rnaseq_pa14_tpm_log <- hlog(rnaseq_pa14_tpm)

# load SRA strain annotations ---------------------------------------------

run_table <- read.csv(snakemake@input[['sraruntab']], stringsAsFactors = F)
#run_table <- read.csv("inputs/original_compendia/SraRunTable.csv", stringsAsFactors = F)

# determine strain type from metadata -------------------------------------

# Since the compendium is composed of multiple strains and filtering criteria 
# are based on characteristic distributions across the compendium, we want to be
# able to build those distributions from samples aligned to their matched target
# reference to avoid alignment-based artifacts influencing the distributions. 
# To do this we will scrape all fields of the SRA run table to look for "PAO1" 
# and "PA14" to provide best-guesses of what strain each sample is, if annotated.

# search for PAO1
seq_strain_ann_pao1 <- rowSums(apply(run_table, 2, FUN = function(x){
  sapply(x, function(y){
    grepl('PAO1',y) | grepl('PA01',y) # since PAO1 is sometimes misspelled as PA01
  })
})) > 0
# search for PA14
seq_strain_ann_pa14 <- rowSums(apply(run_table, 2, FUN = function(x){
  sapply(x, function(y){
    grepl('PA14', y)
  })
})) > 0
# consolidate PAO11 and PA14 annotations
seq_strain_ann <- sapply(c(1:nrow(run_table)), function(x){
  if(seq_strain_ann_pao1[x]){
    'PAO1'
  } else if(seq_strain_ann_pa14[x]){
    'PA14'
  } else{'No annotation'}
})
# name the strain annotations to match compendium for easy indexing
names(seq_strain_ann) <- run_table$Experiment
                               
# re-order strain annotations  to path compendium samples
seq_strain_ann_order <- seq_strain_ann[colnames(rnaseq_pao1_tpm)]

# re-label samples that weren't in the SRA run table
names(seq_strain_ann_order) <- colnames(rnaseq_pao1_tpm)
seq_strain_ann_order <- replace_na(seq_strain_ann_order, "No annotation")

# examine core gene expression ---------------------------------------

# Using genes known to be stable in their expression as a normalization metric 
# is commonly used in qRT analysis and other mRNA quantification techniques.
# It can provide internal controls and account for technical biases. The same 
# concept can hold true when we look at gene expression across compendia. We can
# choose housekeeping genes for **P. aeruginosa** based on a publication and 
# examine their expression in the array compendium.

# The following set of HK genes comes from Alqarni et al 2016, J Microbiol Methods. 

## PAO1
hks_pao1 <- sapply(c('ppiD','rpoD','rpoS','proC','recA','rpsL','rho','oprL',
                     'tipA','nadB'), 
                   function(x) name_to_PAO1(x))
# only include hk genes if they are in our compendium data
hks_pao1 <- hks_pao1[hks_pao1 %in% rownames(rnaseq_pao1_tpm_log)]
## PA14
hks_pa14 <- sapply(c('ppiD','rpoD','rpoS','proC','recA','rpsL','rho','oprL',
                     'tipA','nadB'), 
                   function(x) name_to_PA14(x))
hks_pa14 <- hks_pa14[hks_pa14 %in% rownames(rnaseq_pa14_tpm_log)]


# determine filtering criteria --------------------------------------------

# We will use the distributions of zero counts and median hk gene expression 
# across the compendium to determine lower and upper thresholding criteria to 
# exclude outlying samples. Using the PAO1- and PA14-aligned data, we can tailor 
# our thresholds to avoid misalignment artifacts by basing the thresholds for 
# each strain on samples of that matched strain.

all.equal(colnames(rnaseq_pao1_tpm_log), names(seq_strain_ann_order))

# zero counts lower and upper thresholds
pao1_zeros <- get_zeros(rnaseq_pao1_tpm_log[ , seq_strain_ann_order == 'PAO1'], .1, .9)
pa14_zeros <- get_zeros(rnaseq_pa14_tpm_log[ , seq_strain_ann_order == 'PA14'], .1, .9)
# median hk expression lower and upper thresholds
pao1_hk <- get_hks(rnaseq_pao1_tpm_log[ , seq_strain_ann_order == 'PAO1'], hks_pao1, .2, .98) #.025,.975 #.2, .98
pa14_hk <- get_hks(rnaseq_pa14_tpm_log[ , seq_strain_ann_order == 'PA14'], hks_pa14, .2, .98)


# apply filtering -------------------------------------------------------

# Now that upper and lower thresholds have been determined, we can apply these 
# criteria across the compendium. I apply the zero-based and hk gene-based as 
# well as the PAO1-based and PA14-based criteria separately and then consolidate 
# in order to be able to re-trace which tests each samples passed and failed.

## PAO1
filt_sp_pao1 <- filter_sparsity(rnaseq_pao1_counts, max_zeros=pao1_zeros[2], min_zeros=pao1_zeros[1])
filt_hk_pao1 <- filter_hks(rnaseq_pao1_tpm_log, hks_pao1, hk_min=pao1_hk[1], hk_max=pao1_hk[2])
filt_pao1_samp_i <- (filt_sp_pao1 & filt_hk_pao1)
## PA14
filt_sp_pa14 <- filter_sparsity(rnaseq_pa14_counts, max_zeros=pa14_zeros[2], min_zeros=pa14_zeros[1])
filt_hk_pa14 <- filter_hks(rnaseq_pa14_tpm_log, hks_pa14, hk_min=pa14_hk[1], hk_max=pa14_hk[2])
filt_pa14_samp_i <- (filt_sp_pa14 & filt_hk_pa14)
# consolidate to include samples that pass PAO1 or PA14 based criteria
filt_i <- filt_pao1_samp_i | filt_pa14_samp_i


# perform filtering -------------------------------------------------------

# For the final compendia we will keep samples that passed either both PAO1 
# zero- and hk-based criteria or both PA14 zero- and hk-based criteria. In this 
# way the PAO1 and PA14 compendia will contain the exact same set of samples, 
# just aligned to difference references.

## PAO1
out_rnaseq_pao1_tpm <- rnaseq_pao1_tpm[rownames(rnaseq_pao1_tpm), filt_i]
out_rnaseq_pao1_counts <- rnaseq_pao1_counts[rownames(rnaseq_pao1_tpm), filt_i]
## PA14
out_rnaseq_pa14_tpm <- rnaseq_pa14_tpm[rownames(rnaseq_pa14_tpm), filt_i]
out_rnaseq_pa14_counts <- rnaseq_pa14_counts[rownames(rnaseq_pa14_tpm), filt_i]


# perform normalization ----------------------------------------------------

# After filtering out outlier datasets there is still a wide range of gene 
# expression distributions that vary sample-to-sample. This is in part due to 
# differences in read depth, though not entirely corrected by read depth 
# normalization such as TPM. This problem has been tackled for differential 
# expression analysis and so we employed a normalization methods used in the 
# DESeq2 method of DE analysis. This method, referred to as the median ratios 
# method (MR) calculates per sample size factors based on the ratio of median 
# expression in each sample to the median expression of a pseudo-reference sample
# defined as the geometric mean of all samples in the compendium. 

normalized_counts_pao1 <- MRnorm(out_rnaseq_pao1_counts)
normalized_counts_pa14 <- MRnorm(out_rnaseq_pa14_counts)


# write out data ----------------------------------------------------------

# Save these compendia in TPM, counts snf normalized counts form noting the 
# 'p2' in the title indicating the parameter combination used for filtering. 

#### PAO1
write.csv(out_rnaseq_pao1_tpm, snakemake@output[['tpm_filt_pao1']], quote = F)
write.csv(out_rnaseq_pao1_counts, snakemake@output[['numreads_filt_pao1']], quote = F)
write.csv(normalized_counts_pao1, snakemake@output[['numreads_filt_norm_pao1']], quote = F)
#### PA14
write.csv(out_rnaseq_pa14_tpm, snakemake@output[['tpm_filt_pa14']], quote = F)
write.csv(out_rnaseq_pa14_counts, snakemake@output[['numreads_filt_pa14']], quote = F)
write.csv(normalized_counts_pa14, snakemake@output[['numreads_filt_norm_pa14']], quote = F)
