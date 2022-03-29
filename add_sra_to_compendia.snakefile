import pandas as pd
import os

m = pd.read_csv("inputs/metadata.csv", header = 0)
m = m[m['accession_in_comp'] == False] # filter the metadata to public data that's not in the compendia
m = m[m['pa_in_reads'] == True] # filter the metadata to samples that contained any Pa reads
SRX = list(m['experiment_accession'])
STRAIN = ['pao1', 'pa14']

h = pd.read_csv("inputs/hogan_metadata.csv", header = 0)
SPU = list(h['sample'])

ALL_SAMPLES = SRX + SPU

rule all:
    input:
        expand("outputs/multiqc/logs_{strain}_multiqc_report.html", strain = STRAIN),
        "outputs/filt_norm_compendia/pao1_aligned_compendium_p2_filtered_tpm.csv",
        "outputs/filt_norm_compendia/pao1_aligned_compendium_p2_filtered_num_reads.csv",
        "outputs/filt_norm_compendia/pao1_aligned_compendium_p2_filtered_counts_norm.csv",
        "outputs/filt_norm_compendia/pa14_aligned_compendium_p2_filtered_tpm.csv",
        "outputs/filt_norm_compendia/pa14_aligned_compendium_p2_filtered_num_reads.csv",
        "outputs/filt_norm_compendia/pa14_aligned_compendium_p2_filtered_counts_norm.csv",
        
rule download_pa14_transcriptome:
    output: "inputs/transcriptomes/pa14.cdna.all.fa.gz"
    resources: mem_mb = 1000
    shell:'''
    wget -O {output} ftp://ftp.ensemblgenomes.org/pub/bacteria/release-49/fasta/bacteria_16_collection/pseudomonas_aeruginosa_ucbpp_pa14_gca_000014625/cdna/Pseudomonas_aeruginosa_ucbpp_pa14_gca_000014625.ASM1462v1.cdna.all.fa.gz
    '''

rule download_pao1_transcriptome:
    output: "inputs/transcriptomes/pao1.cdna.all.fa.gz"
    resources: mem_mb = 1000
    shell:'''
    wget -O {output} ftp://ftp.ensemblgenomes.org/pub/bacteria/release-49/fasta/bacteria_5_collection/pseudomonas_aeruginosa_pao1_gca_000006765/cdna/Pseudomonas_aeruginosa_pao1_gca_000006765.ASM676v1.cdna.all.fa.gz
    '''

rule index_transcriptome:
    input: "inputs/transcriptomes/{strain}.cdna.all.fa.gz"
    output: "outputs/t_indxs/{strain}_cdna_k15/info.json"
    params: out_dir = lambda wildcards: "outputs/t_indxs/" + wildcards.strain + "_cdna_k15"
    conda: "envs/salmon.yml"
    threads: 1
    resources:
        mem_mb=8000 
    shell:'''
    salmon index -t {input} -i {params.out_dir} -k 15
    '''

rule fasterq_dump_sra:
    output: "inputs/raw/sra/{srx}.fastq.gz"
    threads: 1
    resources:
        mem_mb=8000
    params: 
        tmp_dir = "tmp/",
        out_dir = "inputs/raw/sra/",
        out_pre = lambda wildcards: "inputs/raw/sra/" + wildcards.srx
    run:
        row = m.loc[m['experiment_accession'] == wildcards.srx]
        srr = row['run'].values[0]
        if not os.path.exists(str(output[0])):
            shell("fasterq-dump {srr} -o {params.out_pre}.fastq --concatenate-reads -t {params.tmp_dir} -f")
            if os.path.exists(str(params.out_pre) + ".fastq"):
                shell("gzip {params.out_pre}.fastq")
        # if SRA download fails, download from the ENA. Relies on ftp links in metadata table
        print("moving on to ENA download")
        if not os.path.exists(str(output[0])) and not os.path.exists(str(params.out_pre) + ".fastq"):
            fastqs = row['fastq_ftp'].values[0]
            fastqs = fastqs.split(";")
            # when single end, just download
            if len(fastqs) == 1:
                fastq = fastqs[0]
                shell("wget -O {output} ftp://{fastq}")
            # when paired end, download, interleave, and remove _1 and _2 files.
            else:
                fastq_1 = fastqs[0]
                fastq_2 = fastqs[1]
                if not os.path.exists(params.out_pre + "_1.fastq.gz"):
                    shell("wget -O {params.out_pre}_1.fastq.gz ftp://{fastq_1}")

                if not os.path.exists(params.out_pre + "_2.fastq.gz"):
                    shell("wget -O {params.out_pre}_2.fastq.gz ftp://{fastq_2}")

                shell("reformat.sh in1={params.out_pre}_1.fastq.gz in2={params.out_pre}_2.fastq.gz out={output}")

                if os.path.exists(output):
                    os.remove(params.out_pre + "_1.fastq.gz")
                    os.remove(params.out_pre + "_2.fastq.gz")

rule salmon:
    input:
        idx = "outputs/t_indxs/{strain}_cdna_k15/info.json",
        reads = "inputs/raw/sra/{srx}.fastq.gz"
    output: 
        quant = "outputs/salmon/{strain}/{srx}/quant.sf",
        logm = "outputs/salmon/{strain}/{srx}/aux_info/meta_info.json",
    params: 
        indx_dir = lambda wildcards: "outputs/t_indxs/" + wildcards.strain + "_cdna_k15",
        out_dir  = lambda wildcards: "outputs/salmon/" + wildcards.strain + "/" + wildcards.srx 
    conda: "envs/salmon.yml"
    threads: 1
    resources:
        mem_mb=8000 
    shell:'''
    salmon quant -i {params.indx_dir} --softclip --softclipOverhangs \
        --minScoreFraction 0.65 --fldMean 51 --seqBias -l A    \
        -r {input.reads} -o {params.out_dir}
    '''
    
rule salmon_spu:
    """
    This is the same rue as above, with the only difference being that this
    one runs on the lab samples, while the previous one ran on the public SRX
    samples. The rule has to be repeated because the input files differ -- 
    if I had set up my files so that the hogan lab samples and the SRX samples 
    lived in the same folder, then I could use the ALL_SAMPLES wildcard above
    to run this rule only once. however, since my files live in different 
    directories (and I like this organization scheme), I (think) I have to 
    repeat this rule.
    """
    input:
        idx = "outputs/t_indxs/{strain}_cdna_k15/info.json",
        reads = "inputs/raw/hogan/{spu}.fastq.gz"
    output: 
        quant = "outputs/salmon/{strain}/{spu}/quant.sf",
        logm = "outputs/salmon/{strain}/{spu}/aux_info/meta_info.json",
    params: 
        indx_dir = lambda wildcards: "outputs/t_indxs/" + wildcards.strain + "_cdna_k15",
        out_dir  = lambda wildcards: "outputs/salmon/" + wildcards.strain + "/" + wildcards.spu 
    conda: "envs/salmon.yml"
    threads: 1
    resources:
        mem_mb=8000 
    shell:'''
    salmon quant -i {params.indx_dir} --softclip --softclipOverhangs \
        --minScoreFraction 0.65 --fldMean 51 --seqBias -l A    \
        -r {input.reads} -o {params.out_dir}
    '''

rule multiqc_salmon_aux_files:
    """
    compiles all log information for the salmon runs, including library type, reads mapped, etc.
    does not join this file to the compendia logs, 
    """
    input:
        expand("outputs/salmon/{{strain}}/{srx}/aux_info/meta_info.json", srx = SRX),
        expand("outputs/salmon/{{strain}}/{spu}/aux_info/meta_info.json", spu = SPU)
    output: "outputs/multiqc/logs_{strain}_multiqc_report.html"
    params: 
        indir = lambda wildcards: "outputs/salmon/" + wildcards.strain,
        outdir = "outputs/multiqc"
    conda: "envs/multiqc.yml"
    threads: 1
    resources:
        mem_mb=8000 
    shell:'''
    multiqc -o {params.outdir} -i logs_{wildcards.strain} {params.indir}
    '''

rule download_pao1_annotation_map:
    output: "inputs/transcriptomes/pao1_gene_names.csv"
    threads: 1
    resources:
        mem_mb=800 
    shell:'''
    wget -O {output} https://osf.io/ubx4s/download
    '''
    
rule download_p14_annotation_map:
    output: "inputs/transcriptomes/pa14_gene_names.csv"
    threads: 1
    resources:
        mem_mb=800
    shell:'''
    wget -O {output} https://osf.io/ema5c/download
    '''
 
rule combine_quant_sf_files:
    """
    combines the quant sf files into two csv files, one that records ReadNum and
    one that record TPM. Annotations are also switched from the transcriptome index
    to the gene index based on a pre-existing annotation mapping file.
    """
    input:
        annot_map = "inputs/transcriptomes/{strain}_gene_names.csv",
        quant_srx = expand("outputs/salmon/{{strain}}/{srx}/quant.sf", srx = SRX),
        quant_spu = expand("outputs/salmon/{{strain}}/{spu}/quant.sf", spu = SPU),
    output: 
        numreads="outputs/combined_new/num_reads_{strain}.csv",
        tpm="outputs/combined_new/TPM_{strain}.csv"
    conda: "envs/tidyverse.yml"
    threads: 1
    resources:
        mem_mb=4000
    script: "scripts/snakemake_quant_collect.R"
    

rule download_pao1_raw_tpm_compendia:
    output: "inputs/original_compendia/TPM_pao1_cdna_k15.csv"
    threads: 1
    resources:
        mem_mb=800
    shell:'''
    wget -O {output} https://osf.io/urz94/download
    '''

rule download_pao1_raw_numreads_compendia:
    output: "inputs/original_compendia/num_reads_pao1_cdna_k15.csv"
    threads: 1
    resources:
        mem_mb=800
    shell:'''
    wget -O {output} https://osf.io/67j4g/download
    '''

rule download_pa14_raw_tpm_compendia:
    output: "inputs/original_compendia/TPM_pa14_cdna_k15.csv"
    threads: 1
    resources:
        mem_mb=800
    shell:'''
    wget -O {output} https://osf.io/pqkb2/download
    '''

rule download_pa14_raw_numread_compendia:
    output: "inputs/original_compendia/num_reads_pa14_cdna_k15.csv"
    threads: 1
    resources:
        mem_mb=800
    shell:'''
    wget -O {output} https://osf.io/uvpsa/download
    '''

rule combine_compendia_with_new_samples_raw:
    input:
        numreads_new="outputs/combined_new/num_reads_{strain}.csv",
        tpm_new="outputs/combined_new/TPM_{strain}.csv",
        numreads_og="inputs/original_compendia/num_reads_{strain}_cdna_k15.csv",
        tpm_og="inputs/original_compendia/TPM_{strain}_cdna_k15.csv"
    output:
        numreads="outputs/combined_compendia/num_reads_{strain}.csv",
        tpm="outputs/combined_compendia/TPM_{strain}.csv"
    conda: "envs/tidyverse.yml"
    threads: 1
    resources:
        mem_mb=4000
    script: "scripts/snakemake_combine_compendia.R"


rule download_pao1_orthologs:
    output: "inputs/transcriptomes/Pseudomonas_aeruginosa_pao1_orthologs.csv.gz"
    threads: 1
    resources: mem_mb = 1000
    shell:'''
    wget -O {output} https://pseudomonas.com/downloads/pseudomonas/pgd_r_20_2/Pseudomonas_aeruginosa_PAO1_107/Pseudomonas_aeruginosa_PAO1_107_orthologs.csv.gz
    '''

rule download_pa14_orthologs:
    output: "inputs/transcriptomes/Pseudomonas_aeruginosa_pa14_orthologs.csv.gz"
    threads: 1
    resources: mem_mb = 1000
    shell:'''
    wget -O {output} https://pseudomonas.com/downloads/pseudomonas/pgd_r_20_2/Pseudomonas_aeruginosa_UCBPP-PA14_109/Pseudomonas_aeruginosa_UCBPP-PA14_109_orthologs.csv.gz
    '''

rule download_pa14_annotations:
    output: "inputs/transcriptomes/Pseudomonas_aeruginosa_pa14_annotations.csv.gz"
    threads: 1
    resources: mem_mb = 1000
    shell:'''
    wget -O {output} https://pseudomonas.com/downloads/pseudomonas/pgd_r_20_2/Pseudomonas_aeruginosa_UCBPP-PA14_109/Pseudomonas_aeruginosa_UCBPP-PA14_109.csv.gz
    '''

rule download_pao1_annotations:
    output: "inputs/transcriptomes/Pseudomonas_aeruginosa_pao1_annotations.csv.gz"
    threads: 1
    resources: mem_mb = 1000
    shell:'''
wget -O {output} https://pseudomonas.com/downloads/pseudomonas/pgd_r_20_2/Pseudomonas_aeruginosa_PAO1_107/Pseudomonas_aeruginosa_PAO1_107.csv.gz
    '''

rule download_sra_run_table:
    output: "inputs/original_compendia/SraRunTable.csv"
    threads: 1
    resources: mem_mb = 1000
    shell:'''
    wget -O {output} https://osf.io/pt7am/download
    '''

rule filter_and_normalize_compendia:
    input:
        pa14_ortho= "inputs/transcriptomes/Pseudomonas_aeruginosa_pa14_orthologs.csv.gz",
        pao1_ortho= "inputs/transcriptomes/Pseudomonas_aeruginosa_pao1_orthologs.csv.gz",
        pa14_annot= "inputs/transcriptomes/Pseudomonas_aeruginosa_pa14_annotations.csv.gz",
        pao1_annot= "inputs/transcriptomes/Pseudomonas_aeruginosa_pao1_annotations.csv.gz",
        pa14_genes= "inputs/transcriptomes/pa14_gene_names.csv",
        pao1_genes= "inputs/transcriptomes/pao1_gene_names.csv",
        sraruntab = "inputs/original_compendia/SraRunTable.csv",
        annofuncs = "scripts/snakemake_annotation_functions.R",
        filtfuncs = "scripts/filter_functions.R",
        numreads  = expand("outputs/combined_compendia/num_reads_{strain}.csv", strain = STRAIN),
        tpm       = expand("outputs/combined_compendia/TPM_{strain}.csv", strain = STRAIN)
    output:
        tpm_filt_pao1="outputs/filt_norm_compendia/pao1_aligned_compendium_p2_filtered_tpm.csv",
        numreads_filt_pao1="outputs/filt_norm_compendia/pao1_aligned_compendium_p2_filtered_num_reads.csv",
        numreads_filt_norm_pao1="outputs/filt_norm_compendia/pao1_aligned_compendium_p2_filtered_counts_norm.csv",
        tpm_filt_pa14="outputs/filt_norm_compendia/pa14_aligned_compendium_p2_filtered_tpm.csv",
        numreads_filt_pa14="outputs/filt_norm_compendia/pa14_aligned_compendium_p2_filtered_num_reads.csv",
        numreads_filt_norm_pa14="outputs/filt_norm_compendia/pa14_aligned_compendium_p2_filtered_counts_norm.csv",
    conda: "envs/filt_and_norm.yml"
    threads: 1
    resources:
        mem_mb=4000
    script: "scripts/snakemake_filter_and_normalize_compendia.R"


