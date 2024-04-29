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
        filtered_r1 = config['general']['output_dir'] + "filtered/{sample}_R1_001_filtered.fastq",
        filtered_r2 = config['general']['output_dir'] + "filtered/{sample}_R2_001_filtered.fastq"
    output:
        bam = config['general']['output_dir'] + "star_alignment/processed_{sample}_Aligned.sortedByCoord.out.bam",
        bam_index = config['general']['output_dir'] + "star_alignment/processed_{sample}_Aligned.sortedByCoord.out.bam.bai"
    params:
        star_index = config['star']['index'],
        gtf = config['star']['gtf'],
        mismatches = config['star']['mismatches'],
        maxmultiplemapping = config['star']['maxmultiplemapping'],
        runThreadN = 16,
        outFileNamePrefix = lambda wildcards: config['general']['output_dir'] + "star_alignment/processed_" + wildcards.sample + "_"
    threads: 16
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
            --limitOutSJcollapsed 100000 \
            --limitIObufferSize 3000000 3000000 \
            --outFilterType BySJout \
            --alignSJoverhangMin 8 \
            --outTmpKeep None

        samtools index -@ {params.runThreadN} {output.bam}
        """
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
