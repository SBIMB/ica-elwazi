#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process results_projects {
    input:
    val semaphor
    path config_file
    val version_no

    output:
    env "census_project_id", emit: census_project
    env "msvcf_project_id", emit: msvcf_project

    script:
    """
    project_name=\$(cat ${config_file} | jq -r .ica_job_project.name)
    prefix=\${project_name%-jobs}

    census_project="\${prefix}-cohort-census"
    census_project_id=\$(icav2 projects list -o table | grep \${census_project} | cut -f1)

    msvcf_project="\${prefix}-msvcf-version-${version_no}"
    msvcf_project_id=\$(icav2 projects list -o table | grep \${msvcf_project} | cut -f1)
    """
}

process cohort_census {
    input:
    val project_id
    val batch_no

    output:
    path "data/batch-${batch_no}/", emit: dest
    env "dir", emit: src

    script:
    """
    dir="/data/batch-${batch_no}/"
    
    icav2 projectdata download --project-id ${project_id} \${dir}* .
    """
}

process multi_sample_vcf {
    input:
    val project_id
    val version_no

    output:
    path "data/", emit: dest
    env "dir", emit: src

    script:
    """
    dir="/data/"
    
    icav2 projectdata download --project-id ${project_id} \${dir}* .
    """
}

workflow {

    main:
    // For testing purposes
    results_projects('ok', "/Users/baruch/code/ica-igg-pipeline/test_config.json", 1)

    view(results_projects.out)

    cohort_census(results_projects.out.census_project, 1)
    multi_sample_vcf(results_projects.out.msvcf_project, 1)

    publish:
    consensus_data = cohort_census.out.dest
    msvcf_data = multi_sample_vcf.out.dest
}

output {
    consensus_data {
        path 'cohort-census'
        mode 'copy'
    }
    msvcf_data {
        path 'msvcf'
        mode 'copy'
    }
}
