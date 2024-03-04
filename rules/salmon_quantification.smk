rule run_salmon_quant:
    conda:
        config['envs']['salmon']
    input:
        filtered_r1 = config['general']['fastq_folder'] + "/filtered/{sample}_R1_001_filtered.fastq.gz",
        filtered_r2 = config['general']['fastq_folder'] + "/filtered/{sample}_R2_001_filtered.fastq.gz"
    output:
        quant = config['general']['output_dir'] + '/salmon_quant/{sample}/quant.sf'
    params:
        index = config['salmon']['index_dir']
    threads: 8
    log:
        salmon_log = "logs/salmon/{wildcards.sample}.log"
    shell:
        """
        mkdir -p {output.quant.rsplit('/', 1)[0]}
        salmon quant \
            --libType A \
            --validateMappings \
            --gcBias \
            --quiet \
            --numGibbsSamples 30 \
            --threads {threads} \
            -i {params.index} \
            -1 {input.filtered_r1} \
            -2 {input.filtered_r2} \
            --output {output.quant.rsplit('/', 1)[0]} &>> {log.salmon_log}
        """
