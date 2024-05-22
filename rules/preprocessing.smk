"""
Preprocessing snakemake file:
1.- fastqc
2.- Trimming
3.- remove contaminants
4.- Summary QC
"""
rule preliminar_fastqc:
    conda:
        config['envs']['fastqc']
    input:
        r1 = config['general']['fastq_folder'] + "{sample}_R1_001.fastq.gz",
        r2 = config['general']['fastq_folder'] + "{sample}_R2_001.fastq.gz"
    output:
        r1_report = config['general']['output_dir'] + "fastqc_reports/{sample}_R1_001_fastqc.html",
        r2_report = config['general']['output_dir'] + "fastqc_reports/{sample}_R2_001_fastqc.html"
    params:
        outdir = config['general']['output_dir'] + "fastqc_reports/"
    log:
        "logs/fastqc/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} && 
        fastqc {input.r1} --outdir {params.outdir} > {log} 2>&1 && 
        fastqc {input.r2} --outdir {params.outdir} >> {log} 2>&1
        """

rule trimgalore:
    conda:
        config['envs']['trimgalore']
    input:
        r1 = config['general']['fastq_folder'] + "{sample}_R1_001.fastq.gz",
        r2 = config['general']['fastq_folder'] + "{sample}_R2_001.fastq.gz"
    output:
        trimmed_r1 = config['general']['output_dir'] + "trimmed/{sample}_R1_001_val_1.fq.gz",
        trimmed_r2 = config['general']['output_dir'] + "trimmed/{sample}_R2_001_val_2.fq.gz"
    params:
        cores = config['trimgalore']['cores'],
        min_length = config['trimgalore']['min_length'],
        min_quality = config['trimgalore']['min_quality'],
        clip_start = config['trimgalore'].get('cut_first_n', 1),
        outdir = config['general']['output_dir']
    log:
        "logs/trimgalore/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}trimmed &&
        mkdir -p {params.outdir}trimmed/fastqc_reports &&
        trim_galore --paired {input.r1} {input.r2} \
        --cores {params.cores} \
        --gzip \
        --fastqc \
        --length {params.min_length} \
        --trim-n \
        --quality {params.min_quality} \
        --clip_R1 {params.clip_start} \
        --clip_R2 {params.clip_start} \
        --fastqc_args "--outdir {params.outdir}/trimmed/fastqc_reports/" \
        --output_dir {params.outdir}/trimmed > {log} 2>&1
        """

rule remove_contaminants:
    conda:
        config['envs']['bowtie']
    input:
        trimmed_r1 = config['general']['output_dir'] + "trimmed/{sample}_R1_001_val_1.fq.gz",
        trimmed_r2 = config['general']['output_dir'] + "trimmed/{sample}_R2_001_val_2.fq.gz"
    output:
        filtered_r1 = config['general']['output_dir'] + "filtered/{sample}_R1_001_filtered.fastq.gz",
        filtered_r2 = config['general']['output_dir'] + "filtered/{sample}_R2_001_filtered.fastq.gz",
        contaminants = config['general']['output_dir'] + "contaminants/{sample}_contaminants.sam"
    params:
        bowtie2_index = config['bowtie2']['index'],
        bowtie2_seedlen = config['bowtie2']['seedlen'],
        outdir = config['general']['output_dir']
    log:
        "logs/remove_contaminants/{sample}.log"
    threads: 6
    shell:
        """
        mkdir -p {params.outdir}filtered &&
        mkdir -p {params.outdir}contaminants &&
        bowtie2 --seedlen={params.bowtie2_seedlen} \
            --threads {threads} \
            --un-conc-gz {params.outdir}filtered/{wildcards.sample} \
            -x {params.bowtie2_index} \
            -1 {input.trimmed_r1} \
            -2 {input.trimmed_r2} \
            -S {output.contaminants} &&
        mv {params.outdir}filtered/{wildcards.sample}.1 {output.filtered_r1} &&
        mv {params.outdir}filtered/{wildcards.sample}.2 {output.filtered_r2} 
        """

rule postfiltering_fastqc:
    conda:
        config['envs']['fastqc']
    input:
        filtered_r1 = config['general']['output_dir'] + "filtered/{sample}_R1_001_filtered.fastq.gz",
        filtered_r2 = config['general']['output_dir'] + "filtered/{sample}_R2_001_filtered.fastq.gz"
    output:
        filtered_r1_report = config['general']['output_dir'] + "filtered/fastqc/{sample}_R1_001_filtered_fastqc.html",
        filtered_r2_report = config['general']['output_dir'] + "filtered/fastqc/{sample}_R2_001_filtered_fastqc.html"
    params:
        outdir = config['general']['output_dir'] + "filtered/fastqc/"
    log:
        "logs/fastqc_postfiltering/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} && 
        fastqc {input.filtered_r1} --outdir {params.outdir} > {log} 2>&1 && 
        fastqc {input.filtered_r2} --outdir {params.outdir} >> {log} 2>&1
        """
