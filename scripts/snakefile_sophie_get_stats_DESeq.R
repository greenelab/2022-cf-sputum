library(DESeq2)

get_DE_stats_DESeq <- function(grp_metadata_file,
                               expression_file,
                               data_type,
                               run_id, 
                               out_file){
  
  # This function performs DE analysis using DESeq.
  # Expression data in expression_file are grouped based on metadata_file
  #
  # Arguments
  # ---------
  # grp_metadata_file: str
  #   File containing mapping between sample id and group
  #
  # expression_file: str
  #   File containing gene expression data
  #
  # data_type: str
  #   Either 'template' or 'simulated' to label saved output file
  #
  # run_id: str
  #   Used as identifier for different simulated experiments
  
  expression_data <- t(read.csv(expression_file, sep="\t", header=TRUE, row.names=1))
  metadata <- read.csv(grp_metadata_file, sep="\t", header=TRUE)
  
  print("Checking sample ordering...")
  print(all.equal(colnames(expression_data), metadata$sample))
  # re order metadata table if it doesn't match to expression_data
  if(!(all.equal(colnames(expression_data), metadata$sample))) {
    print("reordering samples...")
    metadata <- metadata[order(match(metadata$sample, colnames(expression_data))), ]
  }
 
  metadata$group <- as.factor(metadata$group)
  ddset <- DESeqDataSetFromMatrix(expression_data, colData=metadata, design = ~group)
  
  deseq_object <- DESeq(ddset, quiet=TRUE)
  
  # Note parameter settings:
  # `independentFilter=False`: We have turned off the automatic filtering, which
  # filter out those tests from the procedure that have no, or little
  # chance of showing significant evidence, without even looking at their test statistic.
  # Typically, this results in increased detection power at the same experiment-wide
  # type I error, as measured in terms of the false discovery rate.
  # cooksCutoff=True (default): Cook's distance as a diagnostic to tell if a single sample
  # has a count which has a disproportionate impact on the log fold change and p-values.
  # These genes are flagged with an NA in the pvalue and padj columns
  deseq_results <- results(deseq_object, independentFiltering=FALSE)
  
  deseq_results_df <-  as.data.frame(deseq_results)
  
  # Save summary statistics of DEGs
  write.table(deseq_results_df, file = out_file, row.names = T, sep = "\t", quote = F)
}

get_DE_stats_DESeq(grp_metadata_file = unlist(snakemake@input[['grp_metadata_filename']]),
                   expression_file   = unlist(snakemake@input[['data_filename']]),
                   data_type         = unlist(snakemake@params[['data_type']]),
                   run_id            = unlist(snakemake@params[['run_id']]),
                   out_file          = unlist(snakemake@output[['de_stats']]))
