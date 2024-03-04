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
        r1 = config['general']['fastq_folder'] + "/{sample}_R1_001.fastq.gz",
        r2 = config['general']['fastq_folder'] + "/{sample}_R2_001.fastq.gz"
    output:
        r1_report = config['fastqc']['output_dir'] + "/{sample}_R1_001_fastqc.html",
        r2_report = config['fastqc']['output_dir'] + "/{sample}_R2_001_fastqc.html"
    log:
        "logs/fastqc/{sample}.log"
    shell:
        "mkdir -p {config['fastqc']['output_dir']} && "
        "fastqc {input.r1} --outdir {config['fastqc']['output_dir']} > {log} 2>&1 && "
        "fastqc {input.r2} --outdir {config['fastqc']['output_dir']} >> {log} 2>&1"

rule trimgalore:
    conda:
        config['envs']['trimgalore']
    input:
        r1 = config['general']['fastq_folder'] + "/{sample}_R1_001.fastq.gz",
        r2 = config['general']['fastq_folder'] + "/{sample}_R2_001.fastq.gz"
    output:
        trimmed_r1 = config['general']['fastq_folder'] + "/trimmed/{sample}_R1_001.trimmed.fastq.gz",
        trimmed_r2 = config['general']['fastq_folder'] + "/trimmed/{sample}_R2_001.trimmed.fastq.gz"
    params:
        cores = config['trimgalore']['cores'],
        min_length = config['trimgalore']['min_length'],
        clip_start = config['trimgalore'].get('cut_first_n', 0)
    log:
        "logs/trimgalore/{sample}.log"
    shell:
        """
        mkdir -p {config['general']['fastq_folder']}/trimmed &&  
        mkdir -p {config['general']['fastq_folder']}/trimmed/fastqc_reports &&  
        trim_galore --paired {input.r1} {input.r2} 
        --cores {params.cores} 
        --gzip 
        --fastqc 
        --length {params.min_length} 
        --trim-n 
        --clip_R1 {params.clip_start} 
        --clip_R2 {params.clip_start} 
        --fastqc_args \"--outdir {config['general']['fastq_folder']}/trimmed/fastqc_reports/\" 
        --output_dir {config['general']['fastq_folder']}/trimmed > {log} 2>&1
        """
rule remove_contaminants:
    conda:
        config['envs']['bowtie']
    input:
        r1 = config['general']['fastq_folder'] + "/trimmed/{sample}_R1_001.trimmed.fastq.gz",
        r2 = config['general']['fastq_folder'] + "/trimmed/{sample}_R2_001.trimmed.fastq.gz"
    output:
        filtered_r1 = config['general']['fastq_folder'] + "/filtered/{sample}_R1_001_filtered.fastq.gz",
        filtered_r2 = config['general']['fastq_folder'] + "/filtered/{sample}_R2_001_filtered.fastq.gz",
        contaminants = config['general']['fastq_folder'] + "/contaminants/{sample}_contaminants.sam"
    params:
        bowtie2_index = config['bowtie2']['index'],
        bowtie2_seedlen = config['bowtie2']['seedlen']
    threads: 6
    shell:
        """
        mkdir -p {config['general']['fastq_folder']}/filtered && 
        mkdir -p {config['general']['fastq_folder']}/contaminants &&
        bowtie2 --seedlen={params.bowtie2_seedlen} \
            --threads {threads} \
            --un-conc-gz {output.filtered_r1} \
            -x {params.bowtie2_index} \
            -1 {input.r1} \
            -2 {input.r2} \
            -S {output.contaminants}
        """

rule quantification_QC:
    input:
        filtered_r1 = config['general']['fastq_folder'] + "/filtered/{sample}_R1_001_filtered.fastq.gz",
        filtered_r2 = config['general']['fastq_folder'] + "/filtered/{sample}_R2_001_filtered.fastq.gz",
        contaminants = config['general']['fastq_folder'] + "/contaminants/{sample}_contaminants.sam"
    output:
        qc_file = config['general']['output_dir'] + "/qc/{sample}_qc_report.txt"
    shell:
        """
        mkdir -p {config['general']['fastq_folder']}/qc &&

        # Count total reads from both R1 and R2
        tot_reads_r1=$(zcat {input.filtered_r1} | wc -l)
        tot_reads_r2=$(zcat {input.filtered_r2} | wc -l)
        tot_reads=$(( (tot_reads_r1 + tot_reads_r2) / 4 ))
        
        # Initialize QC report
        printf "Contaminant QC run for {wildcards.sample} on %s\n" "$(date)" > {output.qc_file}
        printf "\\t%s\\t%s\\t%s\\n" "READ_TYPE" "READS" "PERCENTAGE" >> {output.qc_file}
        printf "%s\\t%s\\t%s\\t%s\\n" "{wildcards.sample}" "Total" "$tot_reads" "100" >> {output.qc_file}
        
        # Quantify contaminants
        for contaminant_type in tRNA snRNA snoRNA mtDNA rRNA; do
            contaminant_reads=$(samtools view {input.contaminants} | grep -c "$contaminant_type")
            contaminant_percentage=$(awk "BEGIN{{print(($contaminant_reads / $tot_reads) * 100)}}")
            printf '%s\\t%s\\t%s\\t%.2f\\n' "{wildcards.sample}" "$contaminant_type" "$contaminant_reads" "$contaminant_percentage" >> {output.qc_file}
        done

        # Count filtered reads
        filtered_reads=$(zcat {input.filtered_r1} {input.filtered_r2} | wc -l)
        filtered_reads=$((filtered_reads / 4))
        filtered_percentage=$(awk "BEGIN{{print(($filtered_reads / $tot_reads) * 100)}}")
        printf '%s\\t%s\\t%s\\t%.2f\\n' "{wildcards.sample}" "Passed" "$filtered_reads" "$filtered_percentage" >> {output.qc_file}
        """
