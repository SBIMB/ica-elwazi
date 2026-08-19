#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.version = 1

include { project_data as cohort_census } from '../modules/download'
include { project_data as multi_sample_vcf } from '../modules/download'
include { results_projects } from '../modules/download'

workflow {

    main:
    // For testing purposes
    results_projects('ok', file("./test_config.json"), 1)

    view(results_projects.out)

    cohort_census(results_projects.out.census_project, "data/batch-1")
    multi_sample_vcf(results_projects.out.msvcf_project, "data")

    publish:
    consensus_data = cohort_census.out.flatten().collect().map { files -> [label: "version-${params.version}", files: files] }
    msvcf_data = multi_sample_vcf.out.flatten().collect().map { files -> [label: "version-${params.version}", files: files] }
}

output {
    consensus_data {
        path { results ->
            "${results.label}/cohort-census/"
        }
        mode 'copy'
    }
    msvcf_data {
        path { results ->
            "${results.label}/msvcf"
        }
        mode 'copy'
    }
}
