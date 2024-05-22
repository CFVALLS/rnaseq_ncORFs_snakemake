import os

configfile: "config/config.yml"

# Include rules
include: "rules/preprocessing.smk"
include: "rules/star_alignment.smk"
include: "rules/contaminants.smk"
include: "rules/postprocessing.smk"

# Define sample_list
sample_list_path = config['general']['sample_list']
sample_list = [line.strip() for line in open(sample_list_path, 'r')]

rule all:
    input:
         expand([config['general']['output_dir'] + "fastqc_reports/{sample}_R1_001_fastqc.html",
                config['general']['output_dir'] + "fastqc_reports/{sample}_R2_001_fastqc.html",
                # config['general']['output_dir'] + "trimmed/{sample}_R1_001_val_1.fq.gz",
                # config['general']['output_dir'] + "trimmed/{sample}_R2_001_val_2.fq.gz",
                # config['general']['output_dir'] + "filtered/{sample}_R1_001_filtered.fastq.gz",
                # config['general']['output_dir'] + "filtered/{sample}_R2_001_filtered.fastq.gz",
                # config['general']['output_dir'] + "contaminants/{sample}_contaminants.sam",
                # config['general']['output_dir'] + "star_alignment/processed_{sample}_Aligned.sortedByCoord.out.bam",
                # config['general']['output_dir'] + "star_alignment/processed_{sample}_Aligned.sortedByCoord.out.bam.bai",
                # config['general']['output_dir'] + "featureCounts/{sample}_GRCh38.geneCounts.txt",
                # config['general']['output_dir'] + "filtered/fastqc/{sample}_R1_001_filtered_fastqc.html",
                # config['general']['output_dir'] + "filtered/fastqc/{sample}_R2_001_filtered_fastqc.html",
                config['general']['output_dir'] + "salmon/{sample}_quant",
                config['general']['output_dir'] + "multiqc/multiqc_report.html",
                config['general']['output_dir'] + "tximport/SummarizedExperimentObject.RDS",
                ],
        sample=sample_list)


