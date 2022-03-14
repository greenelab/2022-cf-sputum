## ---------------------------
##
## Script name: filter_functions.R
##
## Description:
##  Includes functions defined for filtering gene expression data by zeros and hk genes
##   get_zeros: calculates the values used for thresholds by zeros at percentiles
##   get_hks: calculates the values used for thresholds by hk gene expression at percentiles
##   filter_sparsity: applies a band-pass filter on number of zero count genes per sample
##   filter_hks: applies a band-pass filter on the median hk gene expression per sample
## 
##  It also includes helper functions:
##   hlog: 'safely' tales the log of data replacing -Inf and Nan with 0
##   MRnorm: sample-wise median of ratio normalization on columns
##
##
## Author: Georgia Doing; Taylor Reiter
##
## Date Created: 2021-04-27
## Date Modified: 2022-03-14
##
## Email: Georgia.Doing.GR@Dartmouth.edu
##
## ---------------------------
##
## Notes: Modified to require DESeq2, remove RVenn dependency
## Source: https://raw.githubusercontent.com/georgiadoing/pa-seq-compendia/main/qc_filtering/filter_functions.R  
##
## ---------------------------
library(DESeq2)
source('scripts/fsqn.R')

hlog <- function(df, base=exp(1)){
  # Description: 'safely' tales the log of data replacing -Inf and Nan with 0
  # df: dataframe of gene expression to be logged
  # base: base for log calculation
  # returns: log of df with happy 0s for 0 or neg values
  dfl <- data.frame(log(data.matrix(df),base))
  dfl[dfl == -Inf] <- 0
  dfl[dfl < 0] <- 0
  dfl
}

MRnorm <- function(data){
  # Description: sample-wise median of ratio normalization on columns
  # data: dataframe of gene expression counts to be normalized via DESeq2
  # returns: dataframe normalized to geometric mean of all samples
  meta <- data.frame(experiment = colnames(data))
  dds <- DESeqDataSetFromMatrix(countData = 
                                       ceiling(data), 
                                     colData = meta, design = ~ experiment)
  dds <- estimateSizeFactors(dds)
  counts(dds, normalized=TRUE)
}

get_zeros <- function(data, lt= .1, ut = .9){
  # Description: calculates the values used for thresholds by zeros at percentiles
  # data: a dataframe of gene expression data, counts
  # lt: the percentile used to determine lower bounf
  # ut: the percentile used to determine upper bound
  # returns: vector of lower and upper bounds in that order
  quantile(apply(
    data, 2, 
    FUN = function(x) sum( x  == 0)), probs = c(lt,ut))
}

get_hks <- function(data, hk_genes, lt=.2,ut=.98){
  # Description: calculates the values used for thresholds by median hk gene expression at percentiles
  # lt: the percentile used to determine lower bounf
  # ut: the percentile used to determine upper bound
  # returns: vector of lower and upper bounds in that order
  quantile(apply(data[hk_genes,], 
                 2, FUN = function(x) median(x)), probs = c(lt,ut))
}

filter_hks <- function(data, hk_genes=NA, hk_min=4, hk_max=6){
  # Description: applies a band-pass filter on the median hk gene expression
  # data: dataframe of gene expression with rownames as gene names
  # hk_genes: array of gene names
  # hk_min: float of minimun mean expression, default 0
  # hk_max: float of maximum mean expression, deafult infinity
  # returns: array of bools for filtering columns of data
  if(is.na(hk_genes)){
    #HK genes come from Alqarni et al 2016, J Microbiol Methods. 
    hk_genes <- sapply(c('ppiD','rpoD','rpoS','proC','recA','rpsL','rho','oprL','tipA','nadB','ampC'), function(x) name_to_PAO1(x))
  }
  hk_clean <- hk_genes %in% rownames(data)
  if(sum(hk_clean) < 1){
    print('HK genes not in dataframe rownames')
    break
  } else{
    hk_genes <- hk_genes[hk_clean]
  }
  hk_means <- apply(data[hk_genes,], 2, FUN = function(x) median(x))
  hk_filt <- (hk_means >= hk_min) & (hk_means <= hk_max)
  return(hk_filt)
}

filter_sparsity <- function(data, max_zeros=1000, min_zeros = 200, min_peak=0){
  # Description: applies a band-pass filter on the number of zero count genes a sample can have
  # data: dataframe of gene expression with rownames as gene names
  # max_zeros: int of maximum unexpressed gene per sample, default infinity
  # min_peak: float [0-1] of multiplier for peak cutoff, default 1 (0.15 good)
  # returns: array of bools for filtering columns of data
  if(max_zeros < 0 | min_peak < 0){
    print('numer of zeros and peak coefficinet must be positive float values')
    break
  }
  toosp_test <- colSums(data == 0) < max_zeros
  unsparse_test_z <- colSums(data == 0) > min_zeros
  unsparse_test_p <- sapply(colnames(data), function(x){
    h1 <- hist(as.numeric(data[,x]), plot = F, breaks=50)
    if(h1$density[1] <= min_peak){
      F
    } else {T}
  })
  sparse_filt <- toosp_test & unsparse_test_z #& unsparse_test_p 
  return(sparse_filt)
}

