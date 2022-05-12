import pandas as pd

h = pd.read_csv("inputs/hogan_amplicon_metadata.csv", header = 0)

# some missing 16s samples; define samples separately
SEQ     = ['16s', 'its']
SPU_16S = ['GD105A', 'GD105B', 'GD106', 'GD108',
           'GD120',  'GD122',  'GD124', 'GD145',
           'GD153',  'GD201',  'GD203', 'GD204A',
           'GD220',  'GD233',  'GD239', 'GD239', 'GD243']
SPU_ITS = ['GD105A', 'GD105B', 'GD106', 'GD108',
           'GD120',  'GD122',  'GD124', 'GD125',
           'GD133',  'GD138',  'GD145', #'GD153',  
           'GD201',  'GD203', 'GD204A', 'GD204B',
           'GD220',  'GD233',  'GD239', 'GD239', 'GD243']
# ITS GD153 failed filtering
#
rule all:
    input: expand("outputs/amplicon/{seq}/dada2/summary_table.tsv", seq = SEQ)

rule fastp_adapter_trim:
    """
    requires that input 16s data is located in inputs/raw/hogan_16s
    https://doi.org/10.1128/spectrum.01915-21 recommends light trimming
    as lower quality is associated with genomic GC content, leading
    to skewed taxonomic recall when using stringent trimming.
    """
    output: 
        r1 = 'outputs/amplicon/{seq}/fastp/{spu}_R1.fastp.fq.gz',
        r2 = 'outputs/amplicon/{seq}/fastp/{spu}_R2.fastp.fq.gz',
        json = 'outputs/amplicon/{seq}/fastp/{spu}.json',
        html = 'outputs/amplicon/{seq}/fastp/{spu}.html'
    #conda: 'envs/fastp.yml'
    threads: 1
    params: indir = lambda wildcards: "inputs/raw/hogan_" + wildcards.seq + "/"
    benchmark: "benchmarks/amplicon/{seq}/fastp__{spu}.txt"
    run:
        r1_col = "r1_" + wildcards.seq
        r2_col = "r2_" + wildcards.seq
        row = h.loc[h['sample'] == wildcards.spu]
        r1 = row[r1_col].values[0]
        r2 = row[r2_col].values[0]
        shell("fastp -i {params.indir}{r1} -I {params.indir}{r2} -o {output.r1} -O {output.r2} -q 10 -j {output.json} -h {output.html} -l 50 -c -w {threads}")

#rule dummy_expand_seq:
#    input: expand("outputs/amplicon/{seq}/fastp/{{spu}}.json", seq = SEQ)
#    output: touch("outputs/amplicon/.{spu}_fastp_done.txt")

rule download_silva_database:
    output: "inputs/silva_nr_v138.1_train_set.fa.gz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1
    '''

rule download_unite_database:
    """
    originally downloaded here: https://plutof.ut.ee/#/doi/10.15156/BIO/1280089
    requires user authentication as UNITE polls for usage.
    I uploaded to osf to enable automatic download
    """
    output: "inputs/sh_general_release_s_10.05.2021.tgz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/bhxyp/download
    '''

rule decompress_unite_database:
    input: "inputs/sh_general_release_s_10.05.2021.tgz"
    output: "inputs/sh_general_release_s_10.05.2021/sh_general_release_dynamic_s_10.05.2021.fasta"
    shell:'''
    tar xf {input} -C inputs/
    '''
# Ran the dada2 rules in two steps because there are two separate sets of SPU samples,
# since there were missing 16s samples. If those end up getting sequenced, the two dada2
# rules could be collapsed.
rule dada2_16s:
    input: 
        r1 = expand('outputs/amplicon/16s/fastp/{spu}_R1.fastp.fq.gz', spu = SPU_16S),
        r2 = expand('outputs/amplicon/16s/fastp/{spu}_R2.fastp.fq.gz', spu = SPU_16S),
        silva = "inputs/silva_nr_v138.1_train_set.fa.gz",
        unite = "inputs/sh_general_release_s_10.05.2021/sh_general_release_dynamic_s_10.05.2021.fasta"
    output:
        #quality_plts = "outputs/16s/dada2/quality_plots.pdf",
        error_plts = "outputs/amplicon/16s/dada2/error_plots.pdf",
        asvs = "outputs/amplicon/16s/dada2/ASVs.fa",
        counts = "outputs/amplicon/16s/dada2/ASV_counts.tsv",
        tax = "outputs/amplicon/16s/dada2/ASV_taxonomy.tsv",
        summary = "outputs/amplicon/16s/dada2/summary_table.tsv"
    params: seq = "16s"
    conda: "envs/dada2.yml"
    script: "scripts/snakemake_amplicon_dada2.R"

rule dada2_its:
    input: 
        r1 = expand('outputs/amplicon/its/fastp/{spu}_R1.fastp.fq.gz', spu = SPU_ITS),
        r2 = expand('outputs/amplicon/its/fastp/{spu}_R2.fastp.fq.gz', spu = SPU_ITS),
        silva = "inputs/silva_nr_v138.1_train_set.fa.gz",
        unite = "inputs/sh_general_release_s_10.05.2021/sh_general_release_dynamic_s_10.05.2021.fasta"
    output:
        #quality_plts = "outputs/its/dada2/quality_plots.pdf",
        error_plts = "outputs/amplicon/its/dada2/error_plots.pdf",
        asvs = "outputs/amplicon/its/dada2/ASVs.fa",
        counts = "outputs/amplicon/its/dada2/ASV_counts.tsv",
        tax = "outputs/amplicon/its/dada2/ASV_taxonomy.tsv",
        summary = "outputs/amplicon/its/dada2/summary_table.tsv"
    params: seq = "its"
    conda: "envs/dada2.yml"
    script: "scripts/snakemake_amplicon_dada2.R"

