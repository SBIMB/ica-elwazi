#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process aggregate_census {
    conda './.conda/popgen'

    input:
    val semaphor
    path config_file
    path static_config, name: "project-data/config/*"
    path batches, name: "project-data/batch-gvcf/*"
    path version_batches, name: "project-data/version-batch/*"
    path genders, name: "project-data/version-gender/*"
    val version_no
    val shard_ids

    output:
    path "project-analyses/*"

    script:
    """
    popgen-cli dragen-igg submit --step-name aggregate-census \
        --input-project-data-folder-path project-data \
        --input-project-config-file-path ${config_file} \
        --output-analysis-json-folder-path project-analyses \
        --version-ids "${version_no}" \
        --shard-ids "${shard_ids}" \
        --analysis-instance-tier economy 
    """
}
