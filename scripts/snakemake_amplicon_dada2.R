library(dada2)

# set variables -----------------------------------------------------------

# make output file names
fastp_file_dir <- paste0("outputs/amplicon/", snakemake@params[['seq']], "/fastp")
filt_file_dir  <- paste0("outputs/amplicon/", snakemake@params[['seq']], "/filtered/")
# file names of all the forward reads, and one with the reverse
forward_reads <- sort(c(list.files(path = fastp_file_dir, pattern="_R1.fastp.fq.gz", full.names = T)))
reverse_reads <- sort(c(list.files(path = fastp_file_dir, pattern="_R2.fastp.fq.gz", full.names = T)))

# make sure forward and reverse reads are matched
stopifnot(all.equal(
  gsub("_R1.fastp.fq.gz", "", forward_reads),
  gsub("_R2.fastp.fq.gz", "", reverse_reads)
))

# make sample vector from file names
samples <- gsub("_R1.fastp.fq.gz", "", basename(forward_reads))

# generate file names for filtered forward and reverse reads
filtered_forward_reads <- paste0(filt_file_dir, samples, "_R1.filtered.fq.gz")
filtered_reverse_reads <- paste0(filt_file_dir, samples, "_R2.filtered.fq.gz")

print("Done generating file and sample names!")

# diagnostic plots and filtering ------------------------------------------
# pdf(snakemake@output[['quality_plts']])
# plotQualityProfile(forward_reads)
# plotQualityProfile(reverse_reads)
# dev.off()

# print("Done making quality plots!")

# For the forward and reverse reads separately -- 
# maxEE: quality filtering threshold; throw the read away if it is likely to have more than 2 erroneous base calls 
# rm.phix: removes any reads that match the PhiX bacteriophage genome 
# multithread: run in parallel if set to TRUE or if a number is given for # cores
# minLen: minimum length reads we want to keep after trimming
# truncQ: (default) 2, trim all bases after first quality score of 2 in a read. 
# maxN: (default) remove any read containing Ns. 
# truncLen: minimum size to trim the reads to keep the quality scores ~>30.

filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads, 
                              reverse_reads, filtered_reverse_reads, 
                              maxEE=c(2,2), rm.phix=TRUE, 
                              multithread=TRUE, minLen=150, 
                              truncLen=c(220, 180))

print("done filtering reads!")
# generate error model of data --------------------------------------------

err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE)

# Viz 
pdf(snakemake@output[['error_plts']])
plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)
dev.off()
print("done generating error model!")
# dereplicate -------------------------------------------------------------

derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples
print("done dereplicating!")
# infer ASVs --------------------------------------------------------------

# infer true biological sequences
# try pool = TRUE first
# if too computationally expensive,  pool = "pseudo"
dada_forward <- dada(derep_forward, err = err_forward_reads, 
                     multithread = TRUE, pool = TRUE)
dada_reverse <- dada(derep_reverse, err = err_reverse_reads, 
                     multithread = TRUE, pool = TRUE)

print("done inferring ASVs")
# merge forward and reverse reads -----------------------------------------

merged_amplicons <- mergePairs(dada_forward, derep_forward, 
                               dada_reverse, derep_reverse, 
                               trimOverhang=TRUE, minOverlap=50)
print("done  merging amplicons!")
# generate count table 1 ----------------------------------------------------

seqtab <- makeSequenceTable(merged_amplicons)

# chimera detection -------------------------------------------------------

seqtab_nochim <- removeBimeraDenovo(seqtab, multithread=T, verbose=T) 
print(paste0("Bimera: ", sum(seqtab_nochim)/sum(seqtab)))

# where did we lose counts ------------------------------------------------

getN <- function(x) sum(getUniques(x))

# make a little table
summary_tab <- data.frame(row.names = samples,
                          dada2_input = filtered_out[ , 1], 
                          filtered = filtered_out[ , 2], 
                          dada_f = sapply(dada_forward, getN), 
                          dada_r = sapply(dada_reverse, getN), 
                          merged = sapply(merged_amplicons, getN), 
                          nonchim = rowSums(seqtab_nochim), 
                          final_perc_reads_retained = round(rowSums(seqtab_nochim)/filtered_out[,1]*100, 1))


# summary table:
write.table(summary_tab, snakemake@output[['summary']], sep="\t", quote=F, col.names=NA)

# Assign taxonomy ---------------------------------------------------------

if(snakemake@params[['seq']] == "16s"){
  taxa <- assignTaxonomy(seqtab_nochim, snakemake@input[['silva']], multithread=T, tryRC=T)
} else {
  taxa <- assignTaxonomy(seqtab_nochim, snakemake@input[['unite']], multithread=T, tryRC=T)
}
# Output results ----------------------------------------------------------
# fix headers
asv_seqs <- colnames(seqtab_nochim)
asv_headers <- vector(dim(seqtab_nochim)[2], mode="character")

for (i in 1:dim(seqtab_nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, snakemake@output[['asvs']])

# count table:
asv_tab <- t(seqtab_nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, snakemake@output[['counts']], sep="\t", quote=F, col.names=NA)

# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, snakemake@output[['tax']], sep="\t", quote=F, col.names=NA)
