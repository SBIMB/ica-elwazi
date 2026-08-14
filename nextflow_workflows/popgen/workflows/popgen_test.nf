#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process run_cli {
    conda './.conda/popgen'

    output:
    stdout

    script:
    """
      popgen-cli dragen-igg --help
    """
}

workflow {
    run_cli() | view()
}
