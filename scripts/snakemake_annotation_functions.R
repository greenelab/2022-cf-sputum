## ---------------------------
##
## Script name: annotation_functions.R
##
## Description: functions to facilitate translation of one annotation type to another (transcript, gene, ortholog)
##
## Args:
##
## Author: Georgia Doing, Taylor Reiter
##
## Date Created: 2021-09-28
## Date Modified: 2022-03-10
##
## Email: Georgia.Doing.GR@Dartmouth.edu
##
## ---------------------------

# load in annotation files 
# (all originally sources from www.pseudomonas.com)

# source url occurs above variable assignation
# https://pseudomonas.com/downloads/pseudomonas/pgd_r_20_2/Pseudomonas_aeruginosa_PAO1_107/Pseudomonas_aeruginosa_PAO1_107_orthologs.csv.gz
PAO1_orth <- read.csv("inputs/transcriptomes/Pseudomonas_aeruginosa_pao1_orthologs.csv.gz")
# https://pseudomonas.com/downloads/pseudomonas/pgd_r_20_2/Pseudomonas_aeruginosa_UCBPP-PA14_109/Pseudomonas_aeruginosa_UCBPP-PA14_109_orthologs.csv.gz
PA14_orth <- read.csv('inputs/transcriptomes/Pseudomonas_aeruginosa_pa14_orthologs.csv.gz')
# https://pseudomonas.com/downloads/pseudomonas/pgd_r_20_2/Pseudomonas_aeruginosa_UCBPP-PA14_109/Pseudomonas_aeruginosa_UCBPP-PA14_109.csv.gz
PA14_ann <- read.csv("inputs/transcriptomes/Pseudomonas_aeruginosa_pa14_annotations.csv.gz")
# https://pseudomonas.com/downloads/pseudomonas/pgd_r_20_2/Pseudomonas_aeruginosa_PAO1_107/Pseudomonas_aeruginosa_PAO1_107.csv.gz
PAO1_ann <- read.csv('inputs/transcriptomes/Pseudomonas_aeruginosa_pao1_annotations.csv.gz')
PAO1_cdna <- read.csv('inputs/transcriptomes/pao1_gene_names.csv')
PA14_cdna <- read.csv('inputs/transcriptomes/pa14_gene_names.csv')

# functions

#' Convert from 3 or 4 letter name to PAO1 number syntax
#' 
#' @param x A 3 or 4 letter gene name
#' @return The corresponding PAO1 number
#' @example 
#' name_to_PAO1('dnaA')
name_to_PAO1 <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PAO1_ann$Locus.Tag[PAO1_ann$Gene.Name == x]
  if(identical(out_num, character(0))){
    return(x)
  }
  return(substring(out_num[1],1,6))
}


PAO1_to_name <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PAO1_ann$Name[substring(PAO1_ann$Locus.Tag,1,6) == x]
  if(identical(out_num, character(0))){
    return(x)
  }
  if(out_num[1] == ''){
    return(x)
  }
  if(grepl(' ',out_num[1])){
    return(x)
  }
  return(out_num[1])
}

name_to_PA14 <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PA14_ann$Locus.Tag[PA14_ann$Gene.Name == x]
  if(identical(out_num, character(0))){
    return(x)
  }
  return(substring(out_num[1],1,10))
}

PA14_to_name <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PA14_ann$Name[substring(PA14_ann$Locus.Tag,1,10) == x]
  if(identical(out_num, character(0))){
    return(x)
  }
  if(out_num == ''){
    return(x)
  }
  return(out_num)
}

PA14_to_PAO1 <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PA14_orth$Locus.Tag..Hit.[PA14_orth$Locus.Tag..Query. == x & PA14_orth$Strain..Hit. == 'Pseudomonas aeruginosa PAO1 (Reference)']
  if(identical(out_num, character(0))){
    return(x)
  }
  return(out_num[1])
}

PAO1_to_PA14 <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PAO1_orth$Locus.Tag..Hit.[PAO1_orth$Locus.Tag..Query. == x & PAO1_orth$Strain..Hit. == 'Pseudomonas aeruginosa UCBPP-PA14']
  if(identical(out_num, character(0))){
    return(x)
  }
  return(out_num[1])
} 

PAO1_to_cDNA <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PAO1_cdna$X1[PAO1_cdna$X2 == x]
  if(identical(out_num, character(0))){
    return(x)
  }
  return(out_num[1])
}

PA14_to_cDNA <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PA14_cdna$X1[PA14_cdna$X2 == x]
  if(identical(out_num, character(0))){
    return(x)
  }
  return(out_num[1])
}

cDNA_to_PAO1 <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PAO1_cdna$X2[PAO1_cdna$X1 == x]
  if(identical(out_num, character(0))){
    return(x)
  }
  return(out_num[1])
}

cDNA_to_PA14 <- function(x){
  if(x==''){
    return(x)
  }
  out_num <- PA14_cdna$X2[PA14_cdna$X1 == x]
  if(identical(out_num, character(0))){
    return(x)
  }
  return(out_num[1])
}