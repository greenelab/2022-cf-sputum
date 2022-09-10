# This snakefile relies on the outputs from add_new_to_compendia.snakefile,
# specifically the file outputs/filt_norm_compendia/*compendium_p2_filtered_num_reads.csv
# It uses this file as well as sets of samples defined in the R script to perform ADAGE analysis.

STRAIN = ['pao1', 'pa14'] # only run pao1 since this is what the eADAGE model is based on
CONTROLS = ['asm', 'm63', 'lb']
GROUPS = ['spu', 'spu_m', 'spu_pub']

rule install_adage:
    output: install = "outputs/adage/adage_installed.txt"
    conda: "envs/adagepath.yml"
    script: "scripts/snakemake_install_ADAGEpath.R"

rule run_adage:
    input:
        install = "outputs/adage/adage_installed.txt",
        num_reads = "outputs/filt_norm_compendia/{strain}_aligned_compendium_p2_filtered_num_reads.csv"
    output: 
        data_activity = "outputs/adage/{strain}_{group}_{control}_data_activity.tsv",
        limma_result = "outputs/adage/{strain}_{group}_{control}_limma_result.tsv",
        volcano_plot = "outputs/adage/{strain}_{group}_{control}_volcano_plot.pdf",
        activity_heatmap_plot = "outputs/adage/{strain}_{group}_{control}_activity_heatmap_plot.pdf",
        marginal_activity = "outputs/adage/{strain}_{group}_{control}_marginal_activity.tsv",
        marginal_limma = "outputs/adage/{strain}_{group}_{control}_marginal_time.tsv",
        signature_overlap_plot = "outputs/adage/{strain}_{group}_{control}_signature_overlap_plot.pdf",
        pathway_association_df = "outputs/adage/{strain}_{group}_{control}_pathway_association_df.tsv",
        combined_activation_and_assoc_df = "outputs/adage/{strain}_{group}_{control}_combined_activation_and_assoc_df.tsv",
        uncharacterized_sigs = "outputs/adage/{strain}_{group}_{control}_uncharacterized_sigs.txt",
        unique_active_sigs_annotated_df = "outputs/adage/{strain}_{group}_{control}_unique_active_sigs_annotated_df.tsv",
    conda: "envs/adagepath.yml"
    script: "scripts/snakemake_adage_two_level.R"
 
