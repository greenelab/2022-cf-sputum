# Using *Pseudomonas aeruginosa* transcriptome profiles to discover facets of cystic fibrosis disease processes

1) Identify metabolite candidates that promote Pa growth 
2) Test whether metabolites that inversely correlate with CF lung function induce Pa virulence-related pathways

## Getting started with this repository


This repository uses conda to manage software installations. 
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).

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

## Tools and data

Recently, a compendia of *P. aeruginosa* gene expression were created using publicly available RNA-seq data (see [preprint](https://doi.org/10.1101/2022.01.24.477642) and [data repository](https://osf.io/s9gyu/)).
This repository leverages these compendia to identify transcriptome signatures relevant to CF. 
Given that new RNA-seq sputum samples have been processed since the compendia construction, this repository includes a snakefile that processes new samples by SRA experiment accession and adds them to the compendia.

We intend to use the software tools eADAGE and SOPHIE to discover and rank metabolite-specific pathways by specific relevance to CF, and sourmash gather to identify specific strains present in each sample.
