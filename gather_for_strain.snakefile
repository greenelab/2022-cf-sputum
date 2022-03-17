import pandas as pd
import os

m = pd.read_csv("inputs/metadata.csv", header = 0)
m = m[m['pa_in_reads'] == True] # filter the metadata to samples that had Pa reads 
SRX = list(m['experiment_accession'])
STRAIN = ['pao1', 'pa14']

h = pd.read_csv("inputs/hogan_metadata.csv", header = 0)
SPU = list(h['sample'])


rule all: 
    input:
        expand("outputs/rnaseq_sourmash_gather/{srx}_gtdb_k31.csv", srx = SRX),
        expand("outputs/rnaseq_sourmash_gather_spu/{spu}_gtdb_k31.csv", spu = SPU)
    
###########################################################################
## Download databases 
###########################################################################

rule sourmash_download_human_sig:
    output: "inputs/sourmash_dbs/GCF_000001405.39_GRCh38.p13_rna.sig"
    resources:
        mem_mb = 1000,
    threads: 1
    shell:'''
    wget -O {output} https://osf.io/anj6b/download
    '''

rule sourmash_download_gtdb_database_k31:
    output: "inputs/sourmash_dbs/gtdb-rs202.genomic.k31.zip" 
    resources:
        mem_mb = 1000,
    threads: 1
    shell:'''
    wget -O {output} https://osf.io/94mzh/download
    '''


###########################################################################
## Publicly available files
###########################################################################

rule rnaseq_sample_sourmash_sketch:
    """
    Because some of the samples analyzed here were not already in the compendia,
    some of them were downloaded already, with the local fastq file located in
    inputs/raw/sra/{srx}.fastq.gz. If that file does exists, this rule will
    sketch from that file. If it does not exists, this rule will streaming download
    the SRA accession and pipe it to the sketch command, that way the raw data 
    doesn't have to be written to disk. Because of this rule operating based on
    the fastq.gz files existing, this snakefile should not be run in parallel with
    the add_sra_to_compendia.snakefile that is also in this directory.
    """
    output: "outputs/rnaseq_sourmash_sketch/{srx}.sig"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000 ,
    threads: 1
    params: 
        infile = lambda wildcards: "inputs/raw/sra/" + wildcards.srx + ".fastq.gz"
    run:
        if os.path.exists(params.infile):   
            shell("sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund -o {output} --name {wildcards.srx} {params.infile}")
        else:
            row = m.loc[m['experiment_accession'] == wildcards.srx]
            srr = row['run'].values[0]
            shell("fastq-dump --disable-multithreading --fasta 0 --skip-technical --readids --read-filter pass --dumpbase --split-spot --clip -Z {srr} |"
                  "sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --name {wildcards.srx} -o {output} -")
    
rule rnaseq_sample_sourmash_gather_against_gtdb:
    input:
        sig="outputs/rnaseq_sourmash_sketch/{srx}.sig",
        db="inputs/sourmash_dbs/gtdb-rs202.genomic.k31.zip",
        human="inputs/sourmash_dbs/GCF_000001405.39_GRCh38.p13_rna.sig"
    output: "outputs/rnaseq_sourmash_gather/{srx}_gtdb_k31.csv"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000 ,
    threads: 1
    benchmark: "benchmarks/rnaseq/sourmash_gather_k31_{srx}.txt"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash gather -o {output} --scaled 2000 -k 31 {input.sig} {input.db} {input.human}
    '''

######################################################################
## Hogan lab samples
######################################################################

rule rnaseq_sourmash_sketch_spu:
    """
    assumes input fastq files are in inputs/raw/hogan
    """
    output: "outputs/rnaseq_sourmash_sketch_spu/{spu}.sig"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000 ,
    threads: 1
    params: indir = "inputs/raw/hogan/"
    benchmark: "benchmarks/rnaseq/sourmash_sketch_{spu}.txt"
    run:
        row = h.loc[h['sample'] == wildcards.spu]
        r1 = row['r1'].values[0]
        r2 = row['r2'].values[0]
        shell("sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund -o {output} --merge {wildcards.spu} {params.indir}{r1} {params.indir}{r2}")


rule rnaseq_sample_sourmash_gather_against_gtdb_spu:
    input:
        sig="outputs/rnaseq_sourmash_sketch_spu/{spu}.sig",
        db="inputs/sourmash_dbs/gtdb-rs202.genomic.k31.zip",
        human="inputs/sourmash_dbs/GCF_000001405.39_GRCh38.p13_rna.sig"
    output: "outputs/rnaseq_sourmash_gather_spu/{spu}_gtdb_k31.csv"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000 ,
    threads: 1
    benchmark: "benchmarks/rnaseq/sourmash_gather_k31_{spu}.txt"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash gather -o {output} --scaled 2000 -k 31 {input.sig} {input.db} {input.human}
    '''
