#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process results_projects {
    errorStrategy 'retry'
    maxRetries 3

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
    icav2 projects list -o table > project_list

    census_project="\${prefix}-cohort-census"
    census_project_id=\$(cat project_list | grep \${census_project} | cut -f1)

    msvcf_project="\${prefix}-msvcf-version-${version_no}"
    msvcf_project_id=\$(cat project_list | grep \${msvcf_project} | cut -f1)
    """
}

process project_data {
    errorStrategy 'retry'
    maxRetries 3

    input:
    val project_id
    val results_path

    output:
    path "*"

    script:
    """
    icav2 projectdata download --project-id ${project_id} /data/* .

    mv ${results_path}/* .
    rm -rf data
    """
}
