# Star Alignment Snakemake file:
# 1.- map reads
# 2.- Strand Detection
# 3.- summerize reads

import os
import glob

rule star_alignment:
    conda:
        config['envs']['star']
    input:
        filtered_r1 = config['general']['output_dir'] + "filtered/{sample}_R1_001_filtered.fastq.gz",
        filtered_r2 = config['general']['output_dir'] + "filtered/{sample}_R2_001_filtered.fastq.gz"
    output:
        bam = config['general']['output_dir'] + "star_alignment/processed_{sample}_Aligned.sortedByCoord.out.bam",
        bam_index = config['general']['output_dir'] + "star_alignment/processed_{sample}_Aligned.sortedByCoord.out.bam.bai"
    params:
        star_index = config['star']['index'],
        gtf = config['star']['gtf'],
        mismatches = config['star']['mismatches'],
        maxmultiplemapping = config['star']['maxmultiplemapping'],
        runThreadN = 10,
        outFileNamePrefix = lambda wildcards: config['general']['output_dir'] + "star_alignment/processed_" + wildcards.sample + "_"
    threads: 10
    shell:
        """
        STAR --genomeDir {params.star_index} \
            --sjdbGTFfile {params.gtf} \
            --runThreadN {params.runThreadN} \
            --runDirPerm All_RWX \
            --twopassMode Basic \
            --readFilesIn {input.filtered_r1} {input.filtered_r2} \
            --readFilesCommand zcat \
            --outFilterMismatchNmax {params.mismatches} \
            --outFilterMultimapNmax {params.maxmultiplemapping} \
            --outSAMattributes All \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            --outFileNamePrefix {params.outFileNamePrefix} \
            --outFilterType BySJout \
            --alignSJoverhangMin 8 \
            --outTmpKeep None

        samtools index -@ {params.runThreadN} {output.bam}
        """

# FeatureCount was not developed to quantify transcript level counts 
rule featureCounts_star:
    input:
        bams= config['general']['output_dir'] + "star_alignment/processed_{sample}_Aligned.sortedByCoord.out.bam"
    output:
        counts= config['general']['output_dir'] + "featureCounts/{sample}_GRCh38.geneCounts.txt"
    conda:
        config['envs']['featurecounts']
    params:
        reference_genome = config['featureCounts']['reference_genome'],
        gtf = config['featureCounts']['gtf'],
        strand_direction = config['featureCounts']['strandness']
    threads: 8
    shell:
        """
        featureCounts -a {params.gtf} \
            -p \
            -o {output.counts} \
            -s {params.strand_direction} \
            -T {threads} \
            -t gene \
            -g gene_id \
            -J \
            -G {params.reference_genome} \
            {input.bams}
        """


# other possibility is to use pseudoalignment of salmon 

# rule salmon_quant:
#     input:
#         bam = config['general']['output_dir'] + "star_alignment/processed_{sample}_Aligned.sortedByCoord.out.bam",
#     output:
#         transcripts_quant = directory(config['general']['output_dir'] + "salmon/{sample}_quant")
#     conda:
#         config['envs']['salmon']
#     params:
#         libType = config['salmon']['strandness'],
#         hh_transcript_fasta = config['salmon']['transcripts']
#     threads: 10
#     shell:"""
#         salmon quant -t {params.hh_transcript_fasta} \
#             -l {params.libType} \
#             -a {input.bam} \
#             -o {output.transcripts_quant} \
#             -p {threads}
#         """

rule salmon_quant_mapping:
    conda:
        config['envs']['salmon']
    input:
        filtered_r1=config['general']['output_dir'] + "filtered/{sample}_R1_001_filtered.fastq.gz",
        filtered_r2=config['general']['output_dir'] + "filtered/{sample}_R2_001_filtered.fastq.gz"
    output:
        transcripts_quant=directory(config['general']['output_dir'] + "salmon/{sample}_quant")
    params:
        libType=config['salmon']['strandness'],
        transcript_index=config['salmon']['index_dir'],
        output_dir=config['general']['output_dir'] + "salmon/",
        mapping_output= config['general']['output_dir'] + "salmon/{sample}.sam"
    threads: 10
    shell:
        """
        salmon quant -i {params.transcript_index} \
                     -l {params.libType} \
                     -1 {input.filtered_r1} -2 {input.filtered_r2} \
                     --validateMappings --seqBias --gcBias \
                     -p {threads} \
                     --writeMappings {params.mapping_output} \
                     -o {output.transcripts_quant}
        """

rule tximport:
    input:
        quant=expand(
           config['general']['output_dir'] + "salmon/{sample}_quant/quant.sf", sample = [line.strip() for line in open(config['general']['sample_list'], 'r')]
        ),
        lib=expand(
            config['general']['output_dir'] + "salmon/{sample}_quant/lib_format_counts.json",
            sample = [line.strip() for line in open(config['general']['sample_list'], 'r')],
        ),
        aux_info=expand(
            config['general']['output_dir'] + "salmon/{sample}_quant/aux_info", sample = [line.strip() for line in open(config['general']['sample_list'], 'r')]
        ),
        cmd_info=expand(
           config['general']['output_dir'] + "salmon/{sample}_quant/cmd_info.json", sample = [line.strip() for line in open(config['general']['sample_list'], 'r')]
        ),
        libparams=expand(
           config['general']['output_dir'] + "salmon/{sample}_quant/libParams", sample = [line.strip() for line in open(config['general']['sample_list'], 'r')]
        ),
        logs=expand(config['general']['output_dir'] + "salmon/{sample}_quant/logs", sample = [line.strip() for line in open(config['general']['sample_list'], 'r')]),
        tx_to_gene= config['tximport']['tx2gene'] # e.g.: "resources/tx2gene.tsv",
    output:
        txi= config['general']['output_dir'] + "tximport/SummarizedExperimentObject.RDS",
    params:
        extra="type='salmon'",
    log:
        "logs/tximport.log"
    wrapper:
        "v3.10.2/bio/tximport"