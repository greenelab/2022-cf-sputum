# import all of the sophie libraries.
# I chose to run all of the python sophie code via the "run" directive, which doesn't allow for the use of an environment.
# I made this choice because it's easier to keep track of everything that is being run, and is more similar to the notebooks that Alex had originally created for running sophie.
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
PUB_COMPARISONS = ['pub-vs-lbpao1', 'pub-vs-lbpa14', 'pub-vs-m63']

# constrain these wildcards so they solve properly
wildcard_constraints:
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
        expand("outputs/sophie/{strain}__{pub_comparison}/generic_gene_summary.tsv", strain = STRAIN, pub_comparison = PUB_COMPARISONS)


rule format_template_experiments_pub:
    """
    This rule will require the outputs from the snakefile add_new_to_compendia.snakefile.
    """
    input:
        metadata="inputs/metadata.csv",
        counts="outputs/combined_compendia/num_reads_{strain}.csv"
    output:
        num_reads="outputs/sophie_template_experiments/{strain}__{pub_comparison}_num_reads.tsv",
        grps="outputs/sophie_template_experiments/{strain}__{pub_comparison}_groups.tsv",
        ponyo="outputs/sophie_template_experiments/{strain}__{pub_comparison}_ponyo.csv"
    conda: "envs/tidyverse.yml"
    threads: 1
    resources: mem_mb=6000
    script: "scripts/snakemake_sophie_format_template_experiment_public.R"

##########################################################
## SOPHIE -- process and train new VAEs
##########################################################
# current thought on strategy is to have a single sophie config file that only
# controls the things that will be the same between all of the experiments (NN arch, 
# num_simulated_runs, latent_dim, epochs, etc.)
# Other params that are specific to an experiment, like file names, are specified as inputs and outputs


rule normalized_template_experiment_data:
    input:
        raw_template_filename = "outputs/sophie_template_experiments/{strain}__{pub_comparison}_num_reads.tsv",
        raw_compendium_filename = "outputs/sophie_training_compendia/{strain}_compendium.tsv",
        scaler_filename = "outputs/sophie/NN_models/{strain}/data/scaler_transform.pickle"
    output:
        mapped_template_filename = "outputs/sophie/{strain}__{pub_comparison}/data/mapped_template_compendium.tsv",
        normalized_template_filename = "outputs/sophie/{strain}__{pub_comparison}/data/normalized_template_compendium.tsv"
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
        normalized_template_filename = "outputs/sophie/{strain}__{pub_comparison}/data/normalized_template_compendium.tsv",
        scaler_filename = "outputs/sophie/NN_models/{strain}/data/scaler_transform.pickle",
        # more than one file is probably needed for the model, but just one is enough for the DAG to build appropriately.
        m1 = expand("outputs/sophie/NN_models/{{strain}}/models/{NN_architecture}/tybalt_2layer_{latent_dim}latent_decoder_model.h5", NN_architecture = NN_ARCHITECTURE, latent_dim = LATENT_DIM)
    output:
        "outputs/sophie/{strain}__{pub_comparison}/pseudo_experiment/selected_simulated_data_x_{run_id}.txt"
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
            local_dir              = "outputs/sophie/" + wildcards.strain + "__" + wildcards.pub_comparison, # changed meaning of local dir so output files go where I want them to, not in the current working directory. Otherwise they would all go in the same "pseudo_experiment" directory. 
            latent_dim             = sophie_params["latent_dim"],
            run                    = wildcards.run_id)

rule process_template_data:
    input:
        config = "config/sophie_hogan_comparisons.tsv",
        raw_template_filename = "outputs/sophie_template_experiments/{strain}__{pub_comparison}_num_reads.tsv",
        grp_metadata_filename = "outputs/sophie_template_experiments/{strain}__{pub_comparison}_groups.tsv",
    output:
        processed_template_filename = "outputs/sophie/{strain}__{pub_comparison}/data/processed_template_compendium.tsv",
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
        grp_metadata_filename = "outputs/sophie_template_experiments/{strain}__{pub_comparison}_groups.tsv",
        simulated_filename =  "outputs/sophie/{strain}__{pub_comparison}/pseudo_experiment/selected_simulated_data_x_{run_id}.txt"
    output:
        out_simulated_filename =  "outputs/sophie/{strain}__{pub_comparison}/pseudo_experiment/selected_simulated_data_x_{run_id}_processed.txt"
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
        grp_metadata_filename = "outputs/sophie_template_experiments/{strain}__{pub_comparison}_groups.tsv",
        # processed_template_filename:
        data_filename = "outputs/sophie/{strain}__{pub_comparison}/data/processed_template_compendium.tsv"
    output:
        de_stats = "outputs/sophie/DE_stats/DE_stats_template_data_{strain}__{pub_comparison}_real.txt"
    params:
        data_type = "template",
        run_id    = "real"
    conda: "envs/diffex.yml"
    script: "scripts/snakefile_sophie_get_stats_DESeq.R"

rule get_DE_stats_simulated_data:
    input:
        grp_metadata_filename = "outputs/sophie_template_experiments/{strain}__{pub_comparison}_groups.tsv",
        # simulated_data_filename:
        data_filename = "outputs/sophie/{strain}__{pub_comparison}/pseudo_experiment/selected_simulated_data_x_{run_id}_processed.txt" 
    output:
        de_stats = "outputs/sophie/DE_stats/DE_stats_simulated_data_{strain}__{pub_comparison}_{run_id}.txt"
    params:
        data_type = "simulated",
        run_id = lambda wildcards: wildcards.run_id # can't be accessed directly in snakemake as wildcard because same script is used above, where there is no wildcard for the run_id value.
    conda: "envs/diffex.yml"
    script: "scripts/snakefile_sophie_get_stats_DESeq.R"

rule rank_genes:
    input:
        config = "config/sophie_hogan_comparisons.tsv",
        template_stats_filename = "outputs/sophie/DE_stats/DE_stats_template_data_{strain}__{pub_comparison}_real.txt",
        simulated_stats_filenames = expand("outputs/sophie/DE_stats/DE_stats_simulated_data_{{strain}}__{{pub_comparison}}_{run_id}.txt", run_id = RUN_IDS)
    output:
        gene_summary_filename = "outputs/sophie/{strain}__{pub_comparison}/generic_gene_summary.tsv"
    run:
        sophie_params = utils.read_config(input.config)
        
        template_DE_stats, simulated_DE_summary_stats = ranking.process_and_rank_genes_pathways(
            template_stats_filename   = input.template_stats_filename,
            local_dir                 = "outputs/sophie/",
            num_simulated_experiments = sophie_params["num_simulated_runs"],
            project_id                = wildcards.strain + "__" + wildcards.pub_comparison,
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
