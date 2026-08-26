#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process project_settings {
    input:
    path config_file

    output:
    env 'key', emit: key
    env 'project_id', emit: id
    env 'project_name', emit: name

    script:
    """
    key=\$(cat ${config_file} | jq -r .ica_api_key)
    project_id=\$(cat ${config_file} | jq -r .ica_job_project.id)
    project_name=\$(cat ${config_file} | jq -r .ica_job_project.name)
    """
}

process upload_sequence_files {
    errorStrategy 'retry'
    maxRetries 3

    input:
    path local_batches, name: "local-batch-gvcf/*"
    val key
    val project_id
    val project_name
    val version_no
    val batch_no

    output:
    path "batch-gvcf/batch-${batch_no}.gvcfs.csv", emit: batch_files
    env 'dest_dir', emit: upload_dir

    script:
    """
    mkdir batch-gvcf
    batch_file=batch-gvcf/batch-${batch_no}.gvcfs.csv
    dest_dir=/sequence_data/version_${version_no}/batch_${batch_no}/

    while IFS= read -r line; do
      filename=\$(basename \$line)
      dest=\${dest_dir}\${filename}
      icav2 projectdata upload \${line} \${dest} --project-id ${project_id} -k ${key}
      icav2 projectdata upload \${line}.tbi \${dest}.tbi --project-id ${project_id} -k ${key}
      echo "ica://${project_name}\${dest}" >> \${batch_file}
    done < local-batch-gvcf/batch-${batch_no}.gvcfs.csv

    icav2 projectdata upload \${batch_file} /meta/batch-gvcf/ --project-id ${project_id} -k ${key} 
    """
}

process upload_files {
    errorStrategy 'retry'
    maxRetries 3

    input:
    val key
    val project_id
    each path(source)
    val dest

    script:
    """
    icav2 projectdata upload ${source} ${dest} --project-id ${project_id} -k ${key}
    """
}

process upload_directory {
    errorStrategy 'retry'
    maxRetries 3

    input:
    val key
    val project_id
    // source ensures that the relevant files are present in the workdir for this process
    each path(source)
    // search_path specifies which files to upload
    val search_path

    script:
    """
    for file in \$(find -L ${search_path} -type f); do
        path=\$dirname(\$file)
        dest="/data/batch\${file#*batch}"
        icav2 projectdata upload \$file \$dest --project-id ${project_id} -k ${key}
    done
    """
}
