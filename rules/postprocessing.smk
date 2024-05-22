rule multiqc:
    conda:
        config['envs']['multiqc']
    input: 
        project_root_dir = config['general']['wd']
    output: 
        html = config['general']['output_dir'] + "multiqc/multiqc_report.html"
    params:
        output_dir = config['general']['output_dir'] + "multiqc/"
    log:
        "logs/multiqc/multiqc.log"
    shell:
        """
        mkdir -p {params.output_dir} &&
        multiqc {input.project_root_dir} -o {params.output_dir} > {log} 2>&1
        """