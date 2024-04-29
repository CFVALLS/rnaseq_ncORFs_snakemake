"""
Contamination analysis :
1.- kraken
2.- Trimming
3.- remove contaminants
4.- Summary QC
"""

rule kraken2:
    input:
        r1 = config['general']['output_dir'] + "trimmed/{sample}_R1_001_val_1.fq.gz",
        r2 = config['general']['output_dir'] + "trimmed/{sample}_R2_001_val_2.fq.gz",
        db = config['kraken2']['db_path']
    output:
        report = config['general']['output_dir'] + "kraken2/{sample}.report"
    params:
        threads = 10
    shell: """
        kraken2 \
        --use-names \
        --threads {params.threads} \
        --db {input.db} \
        --fastq-input \
        --report {output.report} \
        --gzip-compressed \
        --paired \
        {input.r1} \
        {input.r2}
    """
