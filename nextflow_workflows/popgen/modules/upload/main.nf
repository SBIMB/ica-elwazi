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

workflow sequence_files {
    take:
    file_list
    key
    project_id
    project_name
    version_no
    batch_no

    main:
    files_ch = file_list
        .splitText()
        .map { p -> p.trim() }
        .filter { p -> p }
        .map { p -> tuple(file(p), file("${p}.tbi")) }

    dest_dir = "/sequence_data/version_${version_no}/batch_${batch_no}/"
    uploads = upload_sequence_files(files_ch, key, project_id, dest_dir)

    file_list = uploads
        .combine(project_name)
        .map { dest, proj -> "ica://${proj}${dest}".toString() }
        .collectFile(
            name: "batch-${batch_no}.gvcfs.csv",
            newLine: true,
            sort: true,
        )

    emit:
    upload_dir = dest_dir
    batch_files = file_list
}

process upload_sequence_files {
    maxForks 4
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple path(gvcf), path(tbi)
    val key
    val project_id
    val dest_dir

    output:
    env 'dest'

    script:
    """  
    filename=\$(basename ${gvcf})
    dest=${dest_dir}\${filename}
    icav2 projectdata upload ${gvcf} \${dest} --project-id ${project_id} -k ${key}
    icav2 projectdata upload ${tbi} \${dest}.tbi --project-id ${project_id} -k ${key}
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

    output:
    val true, emit: done

    script:
    """
    for file in \$(find -L ${search_path} -type f); do
        path=\$(dirname \$file)
        dest="/data/batch\${path#*batch}/"
        icav2 projectdata upload \$file \$dest --project-id ${project_id} -k ${key}
    done
    """
}
