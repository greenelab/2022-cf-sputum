# import all of the sophie libraries.
# I chose to run all of the python sophie code via the "run" directive, which doesn't allow for the use of an environment.
# Therefore, all of the dependencies for sophie need to be installed in the environment in which this workflow will be executed.

import os
import sys
import glob
import pickle
# for mac, import matplotlib and change backend
import matplotlib
matplotlib.use('TkAgg')
import pandas as pd
import numpy as np
import scipy.stats as ss
from sklearn import preprocessing
from keras.models import load_model
from ponyo import utils, train_vae_modules, simulate_expression_data
from generic_expression_patterns_modules import process, stats, ranking, new_experiment_process

# I have removed the input and output files from the sophie config file.
# Instead, these files are solved through wildcards, and their paths are specified as inputs and outputs for each of the snakemake rules.
STRAIN = ["pa14", "pao1"]
HOGAN_COMPARISONS = ['asm-vs-asm_m', 'spu-vs-spu_m', 'spu-vs-asm', 'spu-vs-m63', 'spu_m-vs-asm_m', 'asm-vs-m63']

# constrain these wildcards so they solve properly
wildcard_constraints:
#    strain="pa..", # pa and then two characters
#    hogan_comparison="...-vs-.*", # three characters, then -vs-, then any number of any characters
    run_id = "\d+" # constrain to digits only

# A few of the config-specified variables also need to be accessible as global variables in the snakefile.
# They're used to name input/output files correctly.
# I read them in below using sophie's config file reader.
# NUM_SIMULATED_RUNS also needed to be specified as a range from 0:num_simulated_runs to properly solve files names, so that is specified below.
sophie_params = utils.read_config("config/sophie_hogan_comparisons.tsv")
NN_ARCHITECTURE    = sophie_params["NN_architecture"]
LATENT_DIM         = sophie_params["latent_dim"]
RUN_IDS            = list(range(0, sophie_params["num_simulated_runs"]))

rule all:
    input: 
        expand("outputs/sophie/{strain}__{hogan_comparison}/generic_gene_summary.tsv", strain = STRAIN, hogan_comparison = HOGAN_COMPARISONS)

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
        num_reads="outputs/sophie_template_experiments/{strain}__{hogan_comparison}_num_reads.tsv",
        grps="outputs/sophie_template_experiments/{strain}__{hogan_comparison}_groups.tsv",
        ponyo="outputs/sophie_template_experiments/{strain}__{hogan_comparison}_ponyo.csv"
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

rule format_template_experiments_clinical_isolates:
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
        raw_compendium_filename="outputs/sophie_training_compendia/{strain}_compendium.tsv"
    output:
        normalized_compendium_filename = "outputs/sophie/NN_models/{strain}/data/normalized_compendium.tsv",
        scaler_filename = "outputs/sophie/NN_models/{strain}/data/scaler_transform.pickle"
    run:
        # set all seeds to get repeatable VAE models
        process.set_all_seeds()

        # Normalize the compendium (0-1)
        process.normalize_compendium(input.raw_compendium_filename, 
                                     output.normalized_compendium_filename, 
                                     output.scaler_filename)

# NOTE: the function train_vae_modules.train_vae() is imported from ponyo.
#       It currently has an internal base_dir parameter that is controlled by the config file but otherwise not exposed to the user. 
#       This breaks my ability to automate it, as the base_dir is one level above the absolute path of the current working directory,
#       which in this case would be outside of the github directory. 
#       I've reproduced the function below, changing the base_dir specification, for now until this is changed in ponyo.
#       original script: https://github.com/greenelab/ponyo/blob/master/ponyo/train_vae_modules.py#L114
#       issue: https://github.com/greenelab/ponyo/issues/43
rule sophie_train_vae:
    input:
        # still needed to control NN architecture
        config="config/sophie_hogan_comparisons.tsv",
        normalized_compendium_filename = "outputs/sophie/NN_models/{strain}/data/normalized_compendium.tsv",
    output:
        m1 = "outputs/sophie/NN_models/{strain}/models/{NN_architecture}/tybalt_2layer_{latent_dim}latent_decoder_model.h5",
        m2 = "outputs/sophie/NN_models/{strain}/models/{NN_architecture}/tybalt_2layer_{latent_dim}latent_decoder_weights.h5",
        m3 = "outputs/sophie/NN_models/{strain}/models/{NN_architecture}/tybalt_2layer_{latent_dim}latent_encoder_model.h5",
        m4 = "outputs/sophie/NN_models/{strain}/models/{NN_architecture}/tybalt_2layer_{latent_dim}latent_encoder_weights.h5",
        hist = "outputs/sophie/NN_models/{strain}/logs/{NN_architecture}/tybalt_2layer_{latent_dim}latent_hist.svg"
    params: 
        base_dir = "outputs/sophie/NN_models",
        dataset_name = lambda wildcards: wildcards.strain
    run:
        # set all seeds to get repeatable VAE models
        process.set_all_seeds()
        
        from ponyo import vae
        
        def train_vae(config_filename, base_dir, dataset_name, input_data_filename):
            # Read in config variables, used to train the VAE
            sophie_params = utils.read_config(config_filename)
    
            # Read data
            normalized_data = pd.read_csv(input_data_filename, header=0, sep="\t", index_col=0)

            # Train (VAE)
            vae.tybalt_2layer_model(
                sophie_params["learning_rate"], 
                sophie_params["batch_size"],
                sophie_params["epochs"],
                sophie_params["kappa"],
                sophie_params["intermediate_dim"],
                sophie_params["latent_dim"],
                sophie_params["epsilon_std"],
                normalized_data,
                base_dir,
                dataset_name,
                sophie_params["NN_architecture"],
                sophie_params["validation_frac"])

        # Train VAE on new compendium data
        train_vae(config_filename = input.config, 
                  base_dir = params.base_dir,
                  dataset_name = params.dataset_name,
                  input_data_filename = input.normalized_compendium_filename) 


rule normalized_template_experiment_data:
    input:
        raw_template_filename = "outputs/sophie_template_experiments/{strain}__{hogan_comparison}_num_reads.tsv",
        raw_compendium_filename = "outputs/sophie_training_compendia/{strain}_compendium.tsv",
        scaler_filename = "outputs/sophie/NN_models/{strain}/data/scaler_transform.pickle"
    output:
        mapped_template_filename = "outputs/sophie/{strain}__{hogan_comparison}/data/mapped_template_compendium.tsv",
        normalized_template_filename = "outputs/sophie/{strain}__{hogan_comparison}/data/normalized_template_compendium.tsv"
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
    multiple times, once for each num_simulated_runs.
    Note because the way this rule is written, and that NN architecture/latent dim isn't propogated through all output file names,
    if multiple NN architectures and latent dims are specified in the config file, the outputs of those different architectures will be overwritten.
    This means that only one NN architecture/latent dim can be specified.
    I think this is a fine design caveat, as only one model should be used for all of the experiments.
    '''
    input:
        config = "config/sophie_hogan_comparisons.tsv",
        normalized_compendium_filename = "outputs/sophie/NN_models/{strain}/data/normalized_compendium.tsv",
        normalized_template_filename = "outputs/sophie/{strain}__{hogan_comparison}/data/normalized_template_compendium.tsv",
        scaler_filename = "outputs/sophie/NN_models/{strain}/data/scaler_transform.pickle",
        # more than one file is probably needed for the model, but just one is enough for the DAG to build appropriately.
        m1 = expand("outputs/sophie/NN_models/{{strain}}/models/{NN_architecture}/tybalt_2layer_{latent_dim}latent_decoder_model.h5", NN_architecture = NN_ARCHITECTURE, latent_dim = LATENT_DIM)
    output:
        "outputs/sophie/{strain}__{hogan_comparison}/pseudo_experiment/selected_simulated_data_x_{run_id}.txt"
    run:
        sophie_params = utils.read_config(input.config)
        vae_model_dir = "outputs/sophie/NN_models/" + wildcards.strain + "/models/" + sophie_params["NN_architecture"]
        # simulate experiment based on template experiment
        normalized_compendium_df = pd.read_csv(input.normalized_compendium_filename, sep="\t", index_col=0, header=0)
        normalized_template_df = pd.read_csv(input.normalized_template_filename, sep="\t", index_col=0, header=0)

        new_experiment_process.embed_shift_template_experiment(
            normalized_data        = normalized_compendium_df,
            template_experiment    = normalized_template_df,
            vae_model_dir          = vae_model_dir,
            selected_experiment_id = "x", # project_id is already recorded in output file path, put this to something random
            scaler_filename        = input.scaler_filename,
            local_dir              = "outputs/sophie/" + wildcards.strain + "__" + wildcards.hogan_comparison, # changed meaning of local dir so output files go where I want them to, not in the current working directory. Otherwise they would all go in the same "pseudo_experiment" directory. 
            latent_dim             = sophie_params["latent_dim"],
            run                    = wildcards.run_id)

rule process_template_data:
    input:
        config = "config/sophie_hogan_comparisons.tsv",
        raw_template_filename = "outputs/sophie_template_experiments/{strain}__{hogan_comparison}_num_reads.tsv",
        grp_metadata_filename = "outputs/sophie_template_experiments/{strain}__{hogan_comparison}_groups.tsv",
    output:
        processed_template_filename = "outputs/sophie/{strain}__{hogan_comparison}/data/processed_template_compendium.tsv",
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
        config = "config/sophie_hogan_comparisons.tsv",
        grp_metadata_filename = "outputs/sophie_template_experiments/{strain}__{hogan_comparison}_groups.tsv",
        simulated_filename =  "outputs/sophie/{strain}__{hogan_comparison}/pseudo_experiment/selected_simulated_data_x_{run_id}.txt"
    output:
        out_simulated_filename =  "outputs/sophie/{strain}__{hogan_comparison}/pseudo_experiment/selected_simulated_data_x_{run_id}_processed.txt"
    run:
        sophie_params = utils.read_config(input.config)
        method = sophie_params['DE_method']
        count_threshold = sophie_params['count_threshold']

        if method == "deseq":
            stats.process_samples_for_DESeq(
                input.simulated_filename,
                grp_metadata_filename     = input.grp_metadata_filename,
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
        grp_metadata_filename = "outputs/sophie_template_experiments/{strain}__{hogan_comparison}_groups.tsv",
        # processed_template_filename:
        data_filename = "outputs/sophie/{strain}__{hogan_comparison}/data/processed_template_compendium.tsv"
    output:
        de_stats = "outputs/sophie/DE_stats/DE_stats_template_data_{strain}__{hogan_comparison}_real.txt"
    params:
        data_type = "template",
        run_id    = "real"
    conda: "envs/diffex.yml"
    script: "scripts/snakemake_sophie_get_stats_DESeq.R"

rule get_DE_stats_simulated_data:
    input:
        grp_metadata_filename = "outputs/sophie_template_experiments/{strain}__{hogan_comparison}_groups.tsv",
        # simulated_data_filename:
        data_filename = "outputs/sophie/{strain}__{hogan_comparison}/pseudo_experiment/selected_simulated_data_x_{run_id}_processed.txt" 
    output:
        de_stats = "outputs/sophie/DE_stats/DE_stats_simulated_data_{strain}__{hogan_comparison}_{run_id}.txt"
    params:
        data_type = "simulated",
        run_id = lambda wildcards: wildcards.run_id # can't be accessed directly in snakemake as wildcard because same script is used above, where there is no wildcard for the run_id value.
    conda: "envs/diffex.yml"
    script: "scripts/snakefile_sophie_get_stats_DESeq.R"

rule rank_genes:
    input:
        config = "config/sophie_hogan_comparisons.tsv",
        template_stats_filename = "outputs/sophie/DE_stats/DE_stats_template_data_{strain}__{hogan_comparison}_real.txt",
        simulated_stats_filenames = expand("outputs/sophie/DE_stats/DE_stats_simulated_data_{{strain}}__{{hogan_comparison}}_{run_id}.txt", run_id = RUN_IDS)
    output:
        gene_summary_filename = "outputs/sophie/{strain}__{hogan_comparison}/generic_gene_summary.tsv"
    run:
        sophie_params = utils.read_config(input.config)
        
        template_DE_stats, simulated_DE_summary_stats = ranking.process_and_rank_genes_pathways(
            template_stats_filename   = input.template_stats_filename,
            local_dir                 = "outputs/sophie/",
            num_simulated_experiments = sophie_params["num_simulated_runs"],
            project_id                = wildcards.strain + "__" + wildcards.hogan_comparison,
            analysis_type             = "DE",
            col_to_rank_by            = sophie_params["rank_genes_by"],
            logFC_name                = sophie_params["DE_logFC_name"], 
            pvalue_name               = sophie_params["DE_pvalue_name"])

        summary_gene_ranks = ranking.generate_summary_table(
            template_stats_filename    = input.template_stats_filename,
            template_ranking_summary  = template_DE_stats,
            simulated_ranking_summary = simulated_DE_summary_stats,
            col_to_rank               = sophie_params["rank_genes_by"],
            local_dir                 = "outputs/sophie/",
            pathway_or_gene           = "gene",
            params                    = sophie_params) # NOTES: the argument params may break the snakefile, as that's a special name in snakemake

        
        summary_gene_ranks.to_csv(output.gene_summary_filename, sep="\t")
