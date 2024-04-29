rule run_salmon_quant:
    conda:
        config['envs']['salmon']
    input:
        filtered_r1 = config['general']['output_dir'] + "/filtered/{sample}_R1_001_filtered.fastq.gz",
        filtered_r2 = config['general']['output_dir'] + "/filtered/{sample}_R2_001_filtered.fastq.gz"
    output:
        quant = config['general']['output_dir'] + '/salmon_quant/{sample}/quant.sf'
    params:
        index = config['salmon']['index_dir'],
        sample_dir = lambda wildcards: config['general']['output_dir'] + '/salmon_quant/' + wildcards.sample,
        strand_direction = config['salmon']['strandness']
    threads: 8
    log:
        salmon_log = "logs/salmon/{sample}.log"
    shell:
        """
        mkdir -p {params.sample_dir}
        salmon quant \
            --libType {params.strand_direction} \
            --validateMappings \
            --gcBias \
            --quiet \
            --numGibbsSamples 30 \
            --threads {threads} \
            -i {params.index} \
            -1 {input.filtered_r1} \
            -2 {input.filtered_r2} \
            -o {params.sample_dir} \
        >> {log.salmon_log} 2>&1
        """
