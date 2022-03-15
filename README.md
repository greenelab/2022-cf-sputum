# Using *Pseudomonas aeruginosa* transcriptome profiles to discover facets of cystic fibrosis (CF) disease processes

1) Identify metabolite candidates that promote Pa growth 
2) Test whether metabolites that inversely correlate with CF lung function induce Pa virulence-related pathways

## Getting started with this repository


This repository uses conda to manage software installations. 
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).
To setup the environment required for this repository, run the following commands:
```
conda env create --name sputum --file environment.yml
conda activate sputum
```


Snakemake can parallelize job submission and modulate resource usage (RAM, CPUs). 
We used the command below in a slurm cluster, but other cluster engines [are also supported](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

```
snakemake -s add_sra_to_compendia.snakefile -j 16 --use-conda --rerun-incomplete --latency-wait 15 -
-resources mem_mb=200000 --cluster "sbatch -t 720 -J comp -p bmm -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -k
```

Alternatively, snakemake can be executed without a cluster:
```
snakemake -s add_sra_to_compendia.snakefile -j 2 --use-conda --rerun-incomplete -k -n
```
These parameters are described below:

+`-s` specifies the snakefile to be executed (here, `add_sra_to_compendia.snakefile`)
+`-j 2` parallelizes the snakefile over two cores (drop to `-j 1` to only run one process at a time)
+ `--use-conda` tells snakemake to use conda to manage software environments for each rule
+ `--rerun-incomplete` tells snakemake to rerun rules when it thinks a file is incomplete, e.g. as may occur if a file is half-finished when a snakemake process is terminated.
+ `-k` indicates for the snakefile to keep running even if a rule fails. Snakemake will attempt to run all rules that don't depend on the output of the failed rule.
+ `-n` specifies a dry run. Remove this to actually execute the snakefile. The dry run is useful to make sure snakemake is running the desired rules for the desired number of times.

## Tools and data

Recently, a compendia of *P. aeruginosa* gene expression were created using publicly available RNA-seq data (see [preprint](https://doi.org/10.1101/2022.01.24.477642), [data repository](https://osf.io/s9gyu/), and [GitHub repository](https://github.com/georgiadoing/pa-seq-compendia)).
This repository leverages these compendia to identify transcriptome signatures relevant to CF. 
Given that new RNA-seq sputum samples have been processed since the construction of the compendia, and since some sputum samples weren't included in the compendia because they were not annotated as containing *P aeruginosa*, this repository includes a snakefile that processes new samples by SRA experiment accession and adds them to the existing compendia.

We intend to use the software tools [eADAGE](https://pubmed.ncbi.nlm.nih.gov/28711280/) and [SOPHIE](https://www.biorxiv.org/content/10.1101/2021.05.24.445440v1) to discover and rank metabolite-specific pathways by specific relevance to CF, and sourmash gather to identify specific strains present in each sample.

## Adding new SRA accessions to the Pa compendia using `add_sra_to_compendia.snakefile`

The `add_sra_to_compendia.snakefile` orchestrates adding new SRX/SRA accessions to the *P aeruginosa* compendia referenced above.
It relies on a metadata sheet, which in this repository is found at `inputs/metadata.csv`.
The pipeline primarily relies on three columns:

1. `run`: SRA run accession (usually begins with ERR, SRR, etc.)
2. `experiment_accession`: SRA experiment accession (usually begins with ERX, SRX, etc)
3. `fastq_ftp`: ftp links for the sequencing data on the European Nucleotide Archive

This file can be created by hand, or a file with all of this information can be downloaded from the European Nucleotide Archive by searching for an SRA run, experiment, or study accession and using the "Download report" function. 

Currently, the pipeline is written to require two additional metadata columns, `accession_in_comp` and `pa_in_reads`, which are used to filter some accessions out before preprocessing them to add to the compendia. I added both columns manually to the metadata table, and each is a logical column designed to subset the metadata quickly. `accession_in_comp` records whether the SRX accession was already in the (filtered and normalized) compendia, and `pa_in_reads` records whether *P. aeruginosa` was detected in the reads at all either by sourmash gather or by salmon quasi mapping).
``` 
m = m[m['accession_in_comp'] == False] # filter the metadata to public data that's not in the compendia
m = m[m['pa_in_reads'] == True] # filter the metadata to samples that contained any Pa reads
```
These columns can be removed from a metadata table and the lines removed from the snakefile and those columns will no longer be required.

Adding new samples to the compendia requires new filtering and normalization, as these processes use the entire data set.
We expect the changes between compendia to be negligible, but will test this soon.
