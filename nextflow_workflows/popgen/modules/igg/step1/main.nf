#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process aggregate_gvcf {
    errorStrategy 'retry'
    maxRetries 3
    conda './.conda/popgen'

    input:
    path config_file
    path static_config, name: "project-data/config/*"
    path batches, name: "project-data/batch-gvcf/*"
    path version_batches, name: "project-data/version-batch/*"
    path genders, name: "project-data/version-gender/*"
    val batch_ids
    val shard_ids

    output:
    path "project-analyses/*"

    script:
    """
    popgen-cli dragen-igg submit --step-name aggregate-gvcf \
        --input-project-data-folder-path project-data \
        --input-project-config-file-path ${config_file} \
        --output-analysis-json-folder-path project-analyses \
        --batch-ids "${batch_ids}" \
        --shard-ids "${shard_ids}" \
        --track-variants false \
        --analysis-instance-tier economy 
    """
}
