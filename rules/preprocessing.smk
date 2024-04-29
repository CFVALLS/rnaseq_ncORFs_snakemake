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
        r1_report = config['fastqc']['output_dir'] + "{sample}_R1_001_fastqc.html",
        r2_report = config['fastqc']['output_dir'] + "{sample}_R2_001_fastqc.html"
    params:
        outdir = config['fastqc']['output_dir']
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
        filtered_r1 = config['general']['output_dir'] + "filtered/{sample}_R1_001_filtered.fastq",
        filtered_r2 = config['general']['output_dir'] + "filtered/{sample}_R2_001_filtered.fastq",
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

# rule quantification_QC:
#     input:
#         filtered_r1 = config['general']['output_dir'] + "filtered/{sample}_R1_001_filtered.fastq",
#         filtered_r2 = config['general']['output_dir'] + "filtered/{sample}_R2_001_filtered.fastq",
#         contaminants = config['general']['output_dir'] + "contaminants/{sample}_contaminants.sam"
#     output:
#         qc_file = config['general']['output_dir'] + "qc/{sample}_qc_report.txt"
#     params:
#         outdir = config['general']['output_dir'] + 'qc'
#     shell:
#         """
#         function count_reads() {{
#             echo $(zcat "$1" | wc -l)
#         }}

#         function calculate_percentage() {{
#             echo $(awk "BEGIN {{print(($1 / $2) * 100)}}")
#         }}

#         mkdir -p {params.outdir}

#         # Count total reads from both R1 and R2
#         tot_reads_r1=$(count_reads "{input.filtered_r1}")
#         tot_reads_r2=$(count_reads "{input.filtered_r2}")
#         tot_reads=$(( (tot_reads_r1 + tot_reads_r2) / 4 ))
        
#         # Initialize QC report
#         echo "Contaminant QC run for {wildcards.sample} on $(date)" > "{output.qc_file}"
#         echo -e "READ_TYPE\tREADS\tPERCENTAGE" >> "{output.qc_file}"
#         echo -e "{wildcards.sample}\tTotal\t$tot_reads\t100" >> "{output.qc_file}"
        
#         # Quantify contaminants
#         for contaminant_type in tRNA snRNA snoRNA mtDNA rRNA; do
#             contaminant_reads=$(samtools view "{input.contaminants}" | grep -c "$contaminant_type")
#             contaminant_percentage=$(calculate_percentage $contaminant_reads $tot_reads)
#             echo -e "{wildcards.sample}\t$contaminant_type\t$contaminant_reads\t$contaminant_percentage" >> "{output.qc_file}"
#         done

#         # Count filtered reads
#         filtered_reads=$(( ( $(count_reads "{input.filtered_r1}") + $(count_reads "{input.filtered_r2}") ) / 4 ))
#         filtered_percentage=$(calculate_percentage $filtered_reads $tot_reads)
#         echo -e "{wildcards.sample}\tPassed\t$filtered_reads\t$filtered_percentage" >> "{output.qc_file}"
#         """
