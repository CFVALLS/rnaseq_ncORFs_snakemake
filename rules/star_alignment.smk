"""
Preprocessing snakemake file:
1.- fastqc
2.- Trimming
3.- remove contaminants
4.- Summary QC
"""
rule star_alignment:
    conda:
        config['envs']['star']
    input:
        filtered_r1 = config['general']['fastq_folder'] + "/filtered/{sample}_R1_001_filtered.fastq.gz",
        filtered_r2 = config['general']['fastq_folder'] + "/filtered/{sample}_R2_001_filtered.fastq.gz"
    output:
        bam = config['general']['output_dir'] + "/star_alignment/processed_{sample}_aligned.sortedByCoord.out.bam",
        bam_index = config['general']['output_dir'] + "/star_alignment/processed_{sample}_aligned.sortedByCoord.out.bam.bai"
    params:
        star_index = config['star']['index'],
        gtf = config['star']['gtf'],
        mismatches = config['star']['mismatches'],
        maxmultiplemapping = config['star']['maxmultiplemapping'],
        runThreadN = 16,
        outFileNamePrefix = lambda wildcards: config['general']['output_dir'] + "/star_alignment/processed_" + wildcards.sample + "_"
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
            --limitOutSJcollapsed 10000000 \
            --limitIObufferSize 300000000 \
            --outFilterType BySJout \
            --alignSJoverhangMin 8 \
            --outTmpKeep None

        samtools index -@ {params.runThreadN} {output.bam}
        """

rule all:
    input:
        counts = config['general']['output_dir'] + "/featureCounts/combined_GRCh38.geneCountsCDS.txt"

rule featureCounts_star:
    conda:
        config['envs']['featurecounts']
    input:
        bams = expand(config['general']['output_dir'] + "/star_alignment/processed_{sample}_aligned.sortedByCoord.out.bam", sample=config['samples'])
    output:
        counts = config['general']['output_dir'] + "/featureCounts/combined_GRCh38.geneCountsCDS.txt"
    params:
        reference_genome = config['featureCounts']['reference_genome'],
        gtf = config['featureCounts']['gtf'],
        strand_direction = config['featureCounts']['strandness']
    threads: 8
    shell:
        """
        featureCounts -a {params.gtf} \
            -o {output.counts} \
            -s {params.strand_direction} \
            -T {threads} \
            -t CDS \
            -g gene_id \
            -J \
            -G {params.reference_genome} \
            {' '.join(input.bams)}
        """
