STRAIN = ["pa14", "pao1"]
HOGAN_COMPARISONS = ['asm-vs-asm_m', 'spu-vs-spu_m', 'spu-vs-asm', 'spu-vs-m63', 'spu_m-vs-asm_m', 'asm-vs-m63']

rule all:
    input: 
        expand("outputs/sophie_training_compendia/{strain}_compendium.tsv", strain = STRAIN),
        expand("outputs/sophie_template_experiments/{strain}_sputum_num_reads.tsv", strain = STRAIN)

###################################################################
## Download existing data products and files
###################################################################

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

#################################################################
## Format files for SOPHIE
#################################################################

rule format_training_compendia:
    input:
        raw_compendium="inputs/original_compendia/num_reads_{strain}_cdna_k15.csv",
        strain="inputs/original_compendia/SRA_annotations.tsv"
    output: 
        raw_compendium="outputs/sophie_training_compendia/{strain}_compendium.tsv",
    conda: "envs/tidyverse.yml"
    threads: 1
    resources: mem_mb=6000
    script: "scripts/snakemake_sophie_format_training_compendium.R"
    

rule format_template_experiments_hogan:
    """
    This rule will require the outputs from the snakefile add_new_to_compendia.snakefile,
    which processes the hogan lab counts into a dataframe. The rule takes this dataframe of
    counts and separates it into six template experiments, described by the
    hogan_comparisons list.
    """
    input:
        metadata="inputs/hogan_metadata.csv",
        counts="outputs/combined_new/num_reads_{strain}.csv"
    output:
        num-reads="outputs/sophie_template_experiments/{strain}-{hogan_comparison}-num-reads.tsv",
        groups="outputs/sophie_template_experiments/{strain}-{hogan_comparison}-groups.tsv",
        ponyo="outputs/sophie_template_experiments/{strain}-{hogan_comparison}-ponyo.csv",
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

##########################################################
## SOPHIE -- process and train new VAEs
##########################################################
# current thought on strategy is to have a single sophie config file that only
# controls the things that will be the same between all of the experiments (NN arch, 
# num_simulated_runs, latent_dim, epochs, etc.)
# Other params that are specific to an experiment, like file names, are specified as inputs and outputs

rule sophie_normalize_compendium:
    input:
        # still needed to control NN architecture
        config="",
        raw_compendium_filename="outputs/sophie_training_compendia/{strain}_compendium.tsv"
    output:
        normalized_compendium_filename = "outputs/sophie/{strain}_{hogan_comparison}/data/normalized_compendium.tsv",
        scaler_filename = "outputs/sophie/{strain}_{hogan_comparison}/data/scaler_transform.pickle",
    params: base_dir = "outputs/sophie"
    run:
        # set all seeds to get repeatable VAE models
        process.set_all_seeds()

        # Normalize the compendium (0-1)
        process.normalize_compendium(input.raw_compendium_filename, 
                                     output.normalized_compendium_filename, 
                                     output.scaler_filename)

rule sophie_train_vae:
    input:
        # still needed to control NN architecture
        config="",
        normalized_compendium_filename = "outputs/sophie/{strain}_{hogan_comparison}/data/normalized_compendium.tsv",
    output:
        m1 = "outputs/sophie/{strain}_{hogan_comparison}/models/{NN_architecture}/tybalt_2layer_{latent_dim}_decoder_model.h5",
        m2 = "outputs/sophie/{strain}_{hogan_comparison}/models/{NN_architecture}/tybalt_2layer_{latent_dim}_decoder_model.h5",
        m3 = "outputs/sophie/{strain}_{hogan_comparison}/models/{NN_architecture}/tybalt_2layer_{latent_dim}_decoder_model.h5",
        m4 = "outputs/sophie/{strain}_{hogan_comparison}/models/{NN_architecture}/tybalt_2layer_{latent_dim}_decoder_model.h5",
        hist = "outputs/sophie/{strain}_{hogan_comparison}/logs/{NN_architecture}/tybalt_2layer_{latent_dim}latent_hist.svg"
    params: base_dir = "outputs/sophie"
    run:
        # set all seeds to get repeatable VAE models
        process.set_all_seeds()

        # Train VAE on new compendium data
        # NEED TO UPDATE THIS WHEN I CAN CONTROL THE BASE_DIR PARAMETER
        train_vae_modules.train_vae(config_filename = input.config, 
                                    input_data_filename = output.normalized_compendium_filename) 

rule normalized_template_experiment_data:
    input:
        raw_template_filename = "outputs/sophie_template_experiments/{strain}_{hogan_comparison}-num-reads.tsv",
        raw_compendium_filename = "outputs/sophie_training_compendia/{strain}_compendium.tsv",
        scaler_filename = "outputs/sophie/{strain}_{hogan_comparison}/data/scaler_transform.pickle"
    output:
        mapped_template_filename = "outputs/sophie/{strain}_{hogan_comparison}/data/mapped_template_compendium.tsv",
        normalized_template_filename = "outputs/sophie/{strain}_{hogan_comparison}/data/normalized_template_compendium.tsv"
    run:
        new_experiment_process.process_template_experiment(
            template_filename = input.raw_template_filename,
            compendium_filename = input.raw_compendium_filename,
            scaler_filename = input.scaler_filename,
            mapped_template_filename = output.mapped_template_filename,
            processed_template_filename = output.normalized_template_filename)


# NOTE: "pseudo_experiment" folder name is hard-coded in new_experiment_process.embed_shift_template_experiment(). 
rule simulate_experiments_based_on_template_experiment:
    '''
    This was originally run as a for loop over num_simulated_runs, but snakemake will coordinate running it 
    multiple times, once for each num_simulated_runs
    '''
    input:
        config = "",
        normalized_compendium_filename = "outputs/sophie/{strain}_{hogan_comparison}/data/normalized_compendium.tsv",
        normalized_template_filename = "outputs/sophie/{strain}_{hogan_comparison}/data/normalized_template_compendium.tsv"
        scaler_filename = "outputs/sophie/{strain}_{hogan_comparison}/data/scaler_transform.pickle"
    output:
        "outputs/sophie/{strain}_{hogan_comparison}/pseudo_experiment/selected_simulated_data_x_{run_id}.txt"
    run:
        sophie_params = utils.read_config(input.config)
        vae_model_dir # PARSE WILDCARDS TO SOLVE, BC IT CONTAINS THE PROJECT ID
        latent_dim    # LATENT DIM IS A ENCODED AS A WILDCARD IN THE WORKFLOW...USE THAT?

        # simulate experiment based on template experiment
        normalized_compendium_df = pd.read_csv(input.normalized_compendium_filename, sep="\t", index_col=0, header=0)
        normalized_template_df = pd.read_csv(input.normalized_template_filename, sep="\t", index_col=0, header=0)

        new_experiment_process.embed_shift_template_experiment(
            normalized_data = normalized_compendium_df,
            template_experiment = normalized_template_df,
            vae_model_dir = vae_model_dir,
            selected_experiment_id = "x", # project_id is already recorded in output file path, put this to something random
            scaler_filename = input.scaler_filename,
            local_dir = "outputs/sophie/" + wildcards.strain + "_" + wildcards.hogan_comparison, # changed meaning of local dir so output files go where I want them to, not in the current working directory. Otherwise they would all go in the same "pseudo_experiment" directory. 
            latent_dim = latent_dim,
            run = wildcards.run_id)


rule process_template_data:
    input:
        raw_template_filename = "outputs/sophie_template_experiments/{strain}_{hogan_comparison}-num-reads.tsv",
        grp_metadata_filename = "outputs/sophie_template_experiments/{strain}-{hogan_comparison}-groups.tsv",
    output:
        processed_template_filename = "outputs/sophie/{strain}_{hogan_comparison}/data/processed_template_compendium.tsv",
    run:
        sophie_params = utils.read_config(input.config)
        method = sophie_params['DE_method']
        count_threshold = sophie_params['count_threshold']

        if method == "deseq":
            stats.process_samples_for_DESeq(
                expression_filename       = input.raw_template_filename,
                grp_metadata_filename     = input.grp_metadata_filename,
                out_expression_filename   = output.processed_template_filename,
                count_threshold           = count_threshold,
                process_metadata_filename = None)


rule process_simulated_data:
    '''
    This was originally run as a for loop over num_simulated_runs, but snakemake will coordinate running it 
    multiple times, once for each num_simulated_runs
    '''
    input:
        grp_metadata_filename = "outputs/sophie_template_experiments/{strain}-{hogan_comparison}-groups.tsv",
        simulated_filename =  "outputs/sophie/{strain}_{hogan_comparison}/pseudo_experiment/selected_simulated_data_x_{run_id}.txt"
    output:
        out_simulated_filename =  "outputs/sophie/{strain}_{hogan_comparison}/pseudo_experiment/selected_simulated_data_x_{run_id}_processed.txt"
    run:
        sophie_params = utils.read_config(input.config)
        method = sophie_params['DE_method']
        count_threshold = sophie_params['count_threshold']

        if method == "deseq":
            stats.process_samples_for_DESeq(
                input.simulated_filename,
                grp_metadata_filename     = input.grp_metadata_filename
                out_expression_filename   = output.out_simulated_filename,
                count_threshold           = count_threshold,
                process_metadata_filename = None)


rule get_DE_stats_real_data:
    '''
    The code for this rule was originally located here: 
    https://github.com/greenelab/generic-expression-patterns/blob/master/generic_expression_patterns_modules/DE_analysis.R#L94
    I've adapted the code to:
    1. run differential expression with the condition as a factor
    2. allow flexible model experiment indication
    3. work within the snakemake automation framework
    Even w/o the changes, the code would need to be reproduced in this repo because it was not installed via pip like the python
    sophie scripts.
    '''
    input:
        grp_metadata_filename = "outputs/sophie_template_experiments/{strain}-{hogan_comparison}-groups.tsv",
        # processed_template_filename:
        data_filename = "outputs/sophie/{strain}_{hogan_comparison}/data/processed_template_compendium.tsv"
    output:
        de_stats = "outputs/sophie/DE_stats/DE_stats_template_data_{strain}_{hogan_comparison}_real.txt"
    params:
        data_type = "template"
        run_id    = "real"
    conda: "envs/diffex.yml"
    script: "scripts/snakefile_sophie_get_stats_DESeq.R"

rule get_DE_stats_simulated_data:
    input:
        grp_metadata_filename = "outputs/sophie_template_experiments/{strain}-{hogan_comparison}-groups.tsv",
        # simulated_data_filename:
        data_filename = "outputs/sophie/{strain}_{hogan_comparison}/pseudo_experiment/selected_simulated_data_x_{run_id}_processed.txt" 
    output:
        de_stats = "outputs/sophie/DE_stats/DE_stats_simulated_data_{strain}_{hogan_comparison}_{run_id}.txt"
    params:
        data_type = "simulated",
        run_id = lambda wildcards: wildcards.run_id # can't be accessed directly in snakemake as wildcard because same script is used above, where there is no wildcard for the run_id value.
    conda: "envs/diffex.yml"
    script: "scripts/snakefile_sophie_get_stats_DESeq.R"

rule rank_genes:
    input:
        config = "",
        template_stats_filename = "outputs/sophie/DE_stats/DE_stats_template_data_{strain}_{hogan_comparison}_real.txt",
        simulated_stats_filenames = expand("outputs/sophie/DE_stats/DE_stats_simulated_data_{{strain}}_{{hogan_comparison}}_{run_id}.txt", run_id = RUN_IDS)
    output:
        gene_summary_filename = "outputs/sophie/{strain}_{hogan_comparison}/generic_gene_summary.tsv"
    run:
        sophie_params = utils.read_config(input.config)

        template_DE_stats, simulated_DE_summary_stats = ranking.process_and_rank_genes_pathways(
            template_stats_filename   = input.template_stats_filename,
            local_dir                 = "outputs/sophie/",
            num_simulated_experiments = num_simulated_runs , # NOT DEFINED ANYWHERE YET
            project_id                = wildcards.strain + "_" + wilcards.hogan_comparison,
            analysis_type             = "DE",
            col_to_rank_by            = "log2FoldChange", # replaces config "rank_genes_by" # THIS IS STILL NEEDED IN THE CONFIG FILE FOR THE NEXT FUNCTION
            logFC_name                = "log2FoldChange", # replaces config "DE_logFC_name"
            pvalue_name               = "padj")           # replaces config "DE_pvalue_name"

        summary_gene_ranks = ranking.generate_summary_table(
            tempate_stats_filename    = input.template_stats_filename,
            template_ranking_summary  = template_DE_stats,
            simulated_ranking_summary = simulated_DE_summary_stats,
            col_to_rank               = "logs2FoldChange", # replaces config "rank_genes_by"
            local_dir                 = "outputs/sophie/",
            pathway_or_gene           = "gene",
            params                    = sophie_params)

        
        summary_gene_ranks.to_csv(output.gene_summary_filename, sep="\t")
