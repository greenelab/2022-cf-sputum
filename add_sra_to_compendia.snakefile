SRX = ['SRX6981205', 'SRX6981204', 'SRX6981203', 'SRX6981202', 'SRX6981201', 
       'SRX6981200', 'SRX6981199', 'SRX6981198', 'SRX6981194', 'SRX6981193', 
       'SRX6981192', 'SRX6981191', 'SRX6981190', 'SRX3789394', 'SRX3789395', 
       'SRX3789396', 'SRX3789397', 'SRX3789390', 'SRX3789391', 'SRX3789392', 
       'ERX2326470', 'ERX2326473', 'ERX2326481', 'ERX2326482', 'SRX7101177', 
       'SRX7101178', 'SRX7101179', 'SRX7101180', 'SRX7101181', 'SRX7101182', 
       'SRX7101183', 'SRX7101184', 'SRX7101185', 'SRX4632310', 'SRX4632311', 
       'SRX4632312', 'SRX11241821', 'SRX11241819', 'SRX11241817']
STRAIN = ['pao1', 'pa14']

rule all:
    input: 
        expand("outputs/salmon/{strain}/{srx}/quant.sf", strain = STRAIN, srx = SRX)
        
rule download_pa14_transcriptome:
    output: "inputs/transcriptomes/pa14.cdna.all.fa.gz"
    shell:'''
    wget -O {output} ftp://ftp.ensemblgenomes.org/pub/bacteria/release-49/fasta/bacteria_16_collection/pseudomonas_aeruginosa_ucbpp_pa14_gca_000014625/cdna/Pseudomonas_aeruginosa_ucbpp_pa14_gca_000014625.ASM1462v1.cdna.all.fa.gz
    '''

rule download_pao1_transcriptome:
    output: "inputs/transcriptomes/pao1.cdna.all.fa.gz"
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
    output: "inputs/raw/sra/{srx}.fastq"
    conda: "envs/sratools.yml"
    threads: 1
    resources:
        mem_mb=8000 
    params: tmp_dir = "tmp/"
    shell:'''
    fasterq-dump {wildcards.srx} -O {output} -t {params.tmp_dir} -f
    '''

rule salmon:
    input:
        idx = "outputs/t_indxs/{strain}_cdna_k15/info.json",
        reads = "inputs/raw/sra/{srx}.fastq"
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
    
#shape_comp/gene_names.R
#shape_comp/quant_collect.R
#shape_comp/logs_collect.py
