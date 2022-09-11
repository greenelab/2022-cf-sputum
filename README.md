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
snakemake -s add_new_to_compendia.snakefile -j 16 --use-conda --rerun-incomplete --latency-wait 15 -
-resources mem_mb=200000 --cluster "sbatch -t 720 -J comp -p bmm -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -k
```

Alternatively, snakemake can be executed without a cluster:
```
snakemake -s add_new_to_compendia.snakefile -j 2 --use-conda --rerun-incomplete -k -n
```
These parameters are described below:

+ `-s` specifies the snakefile to be executed (here, `add_new_to_compendia.snakefile`)
+ `-j 2` parallelizes the snakefile over two cores (drop to `-j 1` to only run one process at a time)
+ `--use-conda` tells snakemake to use conda to manage software environments for each rule
+ `--rerun-incomplete` tells snakemake to rerun rules when it thinks a file is incomplete, e.g. as may occur if a file is half-finished when a snakemake process is terminated.
+ `-k` indicates for the snakefile to keep running even if a rule fails. Snakemake will attempt to run all rules that don't depend on the output of the failed rule.
+ `-n` specifies a dry run. Remove this to actually execute the snakefile. The dry run is useful to make sure snakemake is running the desired rules for the desired number of times.

## Tools and data

Recently, a compendia of *P. aeruginosa* gene expression were created using publicly available RNA-seq data (see [preprint](https://doi.org/10.1101/2022.01.24.477642), [data repository](https://osf.io/s9gyu/), and [GitHub repository](https://github.com/georgiadoing/pa-seq-compendia)).
This repository leverages these compendia to identify transcriptome signatures relevant to CF. 
Given that new RNA-seq sputum samples have been processed since the construction of the compendia, and since some sputum samples weren't included in the compendia because they were not annotated as containing *P aeruginosa*, this repository includes a snakefile that processes new samples by SRA experiment accession and adds them to the existing compendia.

We used the software tools [eADAGE](https://pubmed.ncbi.nlm.nih.gov/28711280/) and [SOPHIE](https://www.biorxiv.org/content/10.1101/2021.05.24.445440v1) to discover and and rank gene expression patterns by specific relevance to CF, and [sourmash gather](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2) to identify specific strains present in each sample.

## Explanation of snakefiles and outputs

Below we include an explanation of the purpose and outputs from each snakefile. 
We uploaded the main outputs from each snakefile to an [Open Science Framework repository](https://osf.io/3mkva).

Variables in the snakefiles:

+ `pa14`: *Pseudomonas aeruginosa* strain PA14. Generally refers to the salmon transcriptome index against which RNA-seq samples was aligned.
+ `pao1`: *Pseudomonas aeruginosa* strain PAO1. Generally refers to the salmon transcriptome index against which RNA-seq samples was aligned.
+ `spu`: *ex vivo* sputum samples from patients with cystic fibrosis that were cultured with *Pseudomonas aeruginosa* in the Hogan lab. 
+ `spu_m`: *ex vivo* sputum samples from patients with cystic fibrosis that had metals added to them and that were cultured with *Pseudomonas aeruginosa* in the Hogan lab.
+ `spu_pub`: Publicly available sputum samples from patients with cystic fibrosis. Samples are detailed in this [metadata file](inputs/metadata.csv). Samples that did not have *P. aeruginosa* detected via `sourmash gather` were labelled as `FALSE` in the column `pa_in_reads` and removed from the analysis.
+ `lb`: LB medium. Publicly available samples that were cultured in LB were used as controls for comparisons against the publicly available sputum samples (`spu_pub`). Selection of which samples to use as controls is detailed in [this notebook](notebooks/20220523_selecting_wt_ctrls.ipynb).
+ `asm`: Artificial sputum medium. *P. aeruginosa* was cultured in the Hogan lab in ASM. RNA-seq samples were used as controls for comparisons against *ex vivo* samples (`spu` and `spu_m`).
+ `asm_m`: Artificial sputum medium treated with metals. *P. aeruginosa* was cultured in the Hogan lab in ASM. RNA-seq samples were used as controls for comparisons against *ex vivo* samples (`spu` and `spu_m`).
+ `m63`: A minimal medium. *P. aeruginosa* was cultured in the Hogan lab in M63. RNA-seq samples were used as controls for comparisons against *ex vivo* samples (`spu` and `spu_m`).

### `add_new_to_compendia.snakefile`: adding new samples to the Pa compendia using

The `add_new_to_compendia.snakefile` orchestrates adding new RNA-seq samples to the *P. aeruginosa* compendia referenced above.
It was used to add new [`SRX` SRA accessions](inputs/metadata.csv) and to add the [RNA-seq samples from the Hogan lab](inputs/hogan_metadata.csv). 
It relies on the metadata sheets `inputs/metadata.csv` and `inputs/hogan_metadata.csv`.

For the new SRA samples, the pipeline primarily relies on three columns:

1. `run`: SRA run accession (usually begins with ERR, SRR, etc.)
2. `experiment_accession`: SRA experiment accession (usually begins with ERX, SRX, etc)
3. `fastq_ftp`: ftp links for the sequencing data on the European Nucleotide Archive

These columns are used to download the publicly available RNA-seq samples from the European Nucleotide Archive and to process them into counts.

The pipeline is written to require two additional metadata columns which are used to filter out accessions that don't contain *P. aeruginosa* (`pa_in_reads`) or that are already in the compendia (`accession_in_comp`). 

``` 
m = m[m['accession_in_comp'] == False] # filter the metadata to public data that's not in the compendia
m = m[m['pa_in_reads'] == True] # filter the metadata to samples that contained any Pa reads
```

For the hogan lab samples, the pipeline primarily relies on the `sample` identifier in the metadata sheet.

Note that adding new samples to the compendia requires new filtering and normalization, as these processes use the entire data set.
You can see how filtering and normalization changed the count values in the compendia in [this notebook](20220315-compare-compendia.ipynb).

This snakefile produces the output folders:

+ `inputs/original_compendia`: downloaded original compendia files
+ `inputs/transcriptomes`: downloaded PA14 and PAO1 cDNA files.
+ `outputs/t_indxs`: salmon index files for PA14 and PAO1
+ `outputs/interleaved_spu`: interleaved Hogan lab RNA-seq samples
+ `outputs/multiqc`: multiqc ran on the salmon `aux_info/meta_info.json` files.
+ `outputs/salmon`: salmon quantified counts for public samples
+ `outputs/salmon_spu`: salmon quantified counts for Hogan lab samples
+ `outputs/combined_compendia`: new combined compendia as transcripts per million and num reads.
+ `outputs/filt_norm_compendia/`: new filtered and normalized compendia.

Note the `_spu` folders contain the Hogan lab samples; because of the input data format varied between public samples and Hogan lab samples, they needed to be separated into different output files to avoid collisions in the snakemake DAG.

### `gather_for_strain.snakefile`: using `sourmash gather` to determine which genomes are detected in the RNA-seq samples

This snakefile relies on outputs from `add_new_to_compendia.snakefile`.
This snakefile determines which genomes are detected in the publicly available sputum samples and in the *ex vivo* cultured Hogan lab samples. 
For more details on `sourmash gather` and to visualize the taxonomic composition of samples, see [this notebook](notebooks/20220316_plot_gather_results.ipynb).

This snakefile produces the output folders:

+ `rnaseq_sourmash_sketch`: sourmash FracMinHash sketches for publicly available sputum samples.
+ `rnaseq_sourmash_sketch_spu`: sourmash FracMinHash sketches for Hogan lab samples.
+ `rnaseq_sourmash_gather`: sourmash gather results for publicly available sputum samples.
+ `rnaseq_sourmashgather_spu`: sourmash gather results for Hogan lab samples.

Note the `_spu` folders contain the Hogan lab samples; because of the input data format varied between public samples and Hogan lab samples, they needed to be separated into different output files to avoid collisions in the snakemake DAG.

### `adage.snakefile`:


### `sophie.snakefile`:
