STRAIN = ["pa14", "pao1"]

rule all:
    input: 
        expand("outputs/sophie_training_compendia/{strain}_compendium.tsv", strain = STRAIN),
        expand("outputs/sophie_template_experiments/{strain}_sputum_num_reads.tsv", strain = STRAIN)

rule download_pao1_raw_numreads_compendia:
    output: "inputs/original_compendia/num_reads_pao1_cdna_k15.csv"
    threads: 1
    resources:
        mem_mb=800
    shell:'''
    wget -O {output} https://osf.io/67j4g/download
    '''

rule download_pa14_raw_numread_compendia:
    output: "inputs/original_compendia/num_reads_pa14_cdna_k15.csv"
    threads: 1
    resources:
        mem_mb=800
    shell:'''
    wget -O {output} https://osf.io/uvpsa/download
    '''

rule download_strain_metadata:
    output: "inputs/original_compendia/SRA_annotations.tsv"
    threads: 1
    resources:
        mem_mb=800
    shell:'''
    wget -O {output} https://raw.githubusercontent.com/greenelab/core-accessory-interactome/master/data/metadata/SRA_annotations.tsv
    '''

rule download_experiment_metadata:
    output: "inputs/original_compendia/SraRunTable.csv"
    threads: 1
    resources: mem_mb = 1000
    shell:'''
    wget -O {output} https://osf.io/pt7am/download
    '''

rule format_training_compendia:
    input:
        compendium="inputs/original_compendia/num_reads_{strain}_cdna_k15.csv",
        strain="inputs/original_compendia/SRA_annotations.tsv"
    output: 
        compendium="outputs/sophie_training_compendia/{strain}_compendium.tsv",
    conda: "envs/tidyverse.yml"
    threads: 1
    resources: mem_mb=6000
    script: "scripts/snakemake_sophie_format_training_compendium.R"
    

rule format_template_experiments_hogan:
    """
    This rule will require the outputs from the snakefile add_new_to_compendia.snakefile,
    which processes the hogan lab counts into a dataframe. The rule takes this dataframe of
    counts and separates it into five template experiments:
        1. sputum samples
        2. sputum samples treated with metals
        3. artificial sputum samples
        4. artificial sputum samples treated with metals
        5. "M" sputum samples (I don't know what these are so i separated them out)
    """
    input:
        compendium="inputs/original_compendia/num_reads_{strain}_cdna_k15.csv",
        metadata="inputs/hogan_metadata.csv",
        counts="outputs/combined_new/num_reads_{strain}.csv"
    output:
        sputum="outputs/sophie_template_experiments/{strain}_sputum_num_reads.tsv",
        meta_sputum = "outputs/sophie_template_experiments/{strain}_sputum_groups.tsv",
        metals="outputs/sophie_template_experiments/{strain}_metals_num_reads.tsv",
        meta_metals = "outputs/sophie_template_experiments/{strain}_metals_groups.tsv",
        ponyo_metals = "outputs/sophie_template_experiments/{strain}_metals_ponyo.csv",
        artificial_sputum ="outputs/sophie_template_experiments/{strain}_sputum_artificial_num_reads.tsv",
        meta_artificial_sputum =  "outputs/sophie_template_experiments/{strain}_sputum_artificial_groups.tsv",
        artificial_metals = "outputs/sophie_template_experiments/{strain}_metals_artificial_num_reads.tsv",
        meta_artificial_metals = "outputs/sophie_template_experiments/{strain}_metals_artificial_groups.tsv",
        m_sputum = "outputs/sophie_template_experiments/{strain}_sputum_m_num_reads.tsv",
        meta_m_sputum = "outputs/sophie_template_experiments/{strain}_sputum_m_groups.tsv",
    conda: "envs/tidyverse.yml"
    threads: 1
    resources: mem_mb=6000
    script: "scripts/snakemake_sophie_format_template_experiment_hogan.R"

rule format_template_experiments_other_sputum:
    """
    This rule will require the outputs from the snakefile add_new_to_compendia.snakefile,
    which processes the hogan lab counts into a dataframe. Use the sputum metadata table 
    to separate counts from SRA experiments (SRX) into different template experiments 
    based on the study they were generated in. Produces 5 total template experiments.
    """

rule formate_template_experiments_clinical_isolates:
    """
    SRA experiments (SRX) that are from clinical isolates are annotated in
    the SRA_annotations.tsv file. Use these annotations, combined with
    BioProject IDs from the SRA Run Table to separate clinical isolate
    samples from different projects into template experiments
    """
