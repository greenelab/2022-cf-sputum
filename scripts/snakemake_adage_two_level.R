library(ADAGEpath)
library(DT)
library(readr)
library(tidyr)
library(plyr)

# This script uses ADAGEpath functions to do comparisons between sets of sputum
# samples and controls to determine differentially activated signatures

# determine control set from snakemake wildcards -------------------------

strain <- snakemake@wildcards[['strain']]
# strain <- "pao1"
controls <- snakemake@wildcards[['control']]
# controls <- "m63"

# make control set against which to do comparisons
if(controls == "asm"){
  control_set <- c('ASM1', 'ASM2', 'ASM3')
} else if(controls == "m63"){
  control_set <- c('M631', 'M632', 'M633')
} else if(controls == "lb" & strain == "pa14"){
  control_set <- c('SRX474156', 'SRX474157', 'SRX8030422', 'SRX8030424',
                   'SRX474130', 'SRX474131', 'SRX7299397', 'SRX7299398',
                   'SRX8030426', 'SRX8030428', 'SRX470383', 'SRX470384',
                   'SRX804118', 'SRX804119', 'SRX4641857', 'SRX4641858')
} else if(controls == "lb" & strain == "pao1"){
  control_set <- c('SRX5661364', 'SRX5661363', 'SRX6976948', 'SRX6976949',
                   'SRX6976950', 'SRX8862509', 'SRX8862510', 'ERX2813653',
                   'ERX2813654', 'ERX2813655', 'ERX2068559', 'ERX2068560',
                   'ERX2068561', 'SRX4579961', 'SRX4579962', 'SRX4579963',
                   'SRX8487076', 'SRX8487077', 'SRX8487078')
} else {
  print("Don't recognize control set")
}


# determine experimental set from snakemake wildcards ---------------------

group <- snakemake@wildcards[['group']]
# group <- "spu_pub"

if(group == "spu"){
  group_set <- c('SPU105', 'SPU106', 'SPU108', 'SPU120', 'SPU122', 'SPU124', 
                 'SPU125', 'SPU133', 'SPU138', 'SPU145', 'SPU153', 'SPU201', 
                 'SPU203', 'SPU204', 'SPU220', 'SPU233', 'SPU239', 'SPU243')
} else if(group == "spu_pub"){
  group_set <- c('SRX6981205', 'SRX6981204', 'SRX6981203', 'SRX6981202', 
                 'SRX6981201', 'SRX6981200', 'SRX6981199', 'SRX6981198', 
                 'SRX6981194', 'SRX6981193', 'SRX6981192', 'SRX6981191',
                 'SRX6981190', 'SRX3789394', 'SRX3789395', 'SRX3789396',
                 'SRX3789397', 'SRX3789390', 'SRX3789391', 'SRX3789392', 
                 'ERX2326469', 'ERX2326470', 'ERX2326471', 'ERX2326472',
                 'ERX2326473', 'ERX2326474', 'ERX2326475', 'ERX2326476', 
                 'ERX2326477', 'ERX2326478', 'ERX2326479', 'ERX2326480', 
                 'ERX2326481', 'ERX2326482', 'ERX2326483', 'SRX7101177', 
                 'SRX7101178', 'SRX7101179', 'SRX7101180', 'SRX7101181', 
                 'SRX7101182', 'SRX7101183', 'SRX7101184', 'SRX7101185', 
                 'SRX4632308', 'SRX4632309', 'SRX4632310', 'SRX4632311', 
                 'SRX4632312')
} else if(group == "spu_m"){
  group_set <- c('SPU105_M', 'SPU106_M', 'SPU108_M', 'SPU120_M', 'SPU122_M', 
                 'SPU124_M', 'SPU125_M', 'SPU133_M', 'SPU138_M', 'SPU145_M', 
                 'SPU153_M', 'SPU201_M', 'SPU203_M', 'SPU204_M', 'SPU220_M', 
                 'SPU233_M', 'SPU239_M', 'SPU243_M')
} else {
  print("comparison group not known")
}

 
# run ADAGE ---------------------------------------------------------------


# Notes
# PAcompendium, eADAGEmodel, and probedistribution are objects that come with the ADAGEpath package.


# Read in normalized and filtered compendium that contains samples of interest.
# Select samples of interest that will be compared with ADAGE.
# The compendium that is read in has been filtered to remove samples with a high
# number of unobserved genes (zeroes), or low expression of house keeping genes.
# The counts are otherwise raw (e.g. produced by salmon), and so are normalized
# with using the ADAGEpath functions
rnaseq_data <- read_csv(snakemake@input[['num_reads']], show_col_types = F) %>%
#rnaseq_data <- read_csv("outputs/filt_norm_compendia/pao1_aligned_compendium_p2_filtered_num_reads.csv", show_col_types = F) %>%
  dplyr::rename(geneID = `...1`) %>%                            # reset index (no colname) to geneID
  dplyr::select(geneID, one_of(control_set), one_of(group_set)) # select control samples, experimental samples

# note that some of the group_set may not be in the RNAseq data because those samples 
# did not pass filtering thresholds when the compendium was originally built

data_raw <- load_dataset(input        = rnaseq_data,
                         isProcessed  = TRUE,
                         isRNAseq     = TRUE,
                         model        = eADAGEmodel,
                         compendium   = PAcompendium,
                         quantile_ref = probedistribution,
                         norm01       = FALSE)

data_normed <- zeroone_norm(input_data = data_raw, use_ref = TRUE, ref_data = PAcompendium)

data_activity <- calculate_activity(input_data = data_normed, model = eADAGEmodel)
write_tsv(data_activity, snakemake@output[['data_activity']])

# build vector for phenotype information. 
# use control set and the final number of columns in the rna_seq data to build this information.
data_pheno <- c(rep("control", length(control_set)), rep(group, (ncol(rnaseq_data) - (length(control_set) + 1))))

limma_result <- build_limma(data_activity,
                            phenotypes = data_pheno,
                            control_pheno = "control",
                            use.bonferroni = TRUE)
write_tsv(limma_result, snakemake@output[['limma_result']])

active_sigs <- get_active_signatures(limma_result = limma_result,
                                     pheno_group  = "both",
                                     method       = "pareto",
                                     N_fronts     = 10) # adjust N fronts to change the number of returned nodes


# plots -------------------------------------------------------------------
pdf(snakemake@output[['volcano_plot']])
plot_volcano(limma_result = limma_result, highlight_signatures = active_sigs,
             interactive = TRUE)
dev.off()


pdf(snakemake@output[['activity_heatmap_plot']])
plot_activity_heatmap(activity = data_activity, signatures = active_sigs)
dev.off()

# overlapping signature removal -------------------------------------------
# check whether active signatures overlap with each other
# plot a heatmap of odds ratios that represent the odds that two sigantures share a specific number of genes
pdf(snakemake@output[['signature_overlap_plot']])
signature_similarity <- plot_signature_overlap(selected_signatures = active_sigs,
                                               model = eADAGEmodel)
dev.off()

# calculate the marginal activities of similar signatures.
# Marginal activity is the activy of sig A after removing genes that it shares with sig B
marginal_activity <- calculate_marginal_activity(input_data = data_normed,
                                                 selected_signatures = active_sigs,
                                                 model = eADAGEmodel)
write_tsv(marginal_activity, snakemake@output[['marginal_activity']])

# build limma model to test whether marginal activities are different between two conditions
marginal_limma <- build_limma(input_data = marginal_activity,
                              phenotypes = data_pheno, control_pheno = "control")
write_tsv(marginal_limma, snakemake@output[['marginal_limma']])


# remove redundant signatures
unique_active_sigs <- remove_redundant_signatures(marginal_limma,
                                                  sig_cutoff = 0.05)
# pathways ----------------------------------------------------------------

KEGG <- fetch_geneset(type = "KEGG")
# we only consider KEGG pathways with more than 3 genes and less than 500 genes
# as meaningful pathways
KEGG_subset <- KEGG[lengths(KEGG) >= 3 & lengths(KEGG) <= 500]

# associate active signatures to KEGG pathways
pathway_association <- annotate_signatures_with_genesets(
  selected_signatures = unique_active_sigs, model = eADAGEmodel,
  genesets = KEGG_subset)
write_tsv(pathway_association, snakemake@output[['pathway_association_df']])

# Calculate the activity of associated pathways inside active signatures
pathway_activity <- signature_geneset_activity(
  signature_geneset_df = pathway_association[, c("signature", "geneset")],
  gene_set_list = KEGG_subset, model = eADAGEmodel, input_data = data_normed)
#plot_activity_heatmap(pathway_activity)

# run a limma test on pathway activities and find pathways that are active
pathway_limma <- build_limma(pathway_activity, phenotypes = data_pheno,
                             control_pheno = "control", use.bonferroni = TRUE)

# combine pathway association and pathway activation test results
combined_result <- combine_geneset_outputs(
  signature_geneset_association = pathway_association,
  geneset_limma_result = pathway_limma)
write_tsv(combined_result, snakemake@output[['combined_activation_and_assoc_df']])

# what signatures are uncharacterized by kegg?
uncharacterized_sigs <- setdiff(unique_active_sigs, pathway_association$signature)
write.table(uncharacterized_sigs, snakemake@output[['uncharacterized_sigs']])

# grab genes in signature and annotate
# loop over nodes in unique_active_sigs, run annotation, and combine into df
unique_active_sigs_annotated_df <- data.frame()
for(signature in unique_active_sigs){
  df <- annotate_genes_in_signatures(selected_signatures = signature,
                                     model = eADAGEmodel,
                                     curated_pathways = KEGG)
  unique_active_sigs_annotated_df <- dplyr::bind_rows(unique_active_sigs_annotated_df, df)
}
write_tsv(unique_active_sigs_annotated_df, snakemake@output[['unique_active_sigs_annotated_df']])



