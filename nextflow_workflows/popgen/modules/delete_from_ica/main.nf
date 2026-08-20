#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process delete_inputs {
    errorStrategy 'retry'
    maxRetries 3

    input:
    val semaphor
    path config_file
    val path

    script:
    """
    project_id=\$(cat ${config_file} | jq -r .ica_job_project.id)

    icav2 projectdata delete ${path} --project-id \${project_id} 
    """
}

process delete_results {
    errorStrategy 'retry'
    maxRetries 3

    input:
    val project_id
    val semaphor

    script:
    """
    icav2 projectdata delete '/data/' --project-id ${project_id}
    """
}
