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
