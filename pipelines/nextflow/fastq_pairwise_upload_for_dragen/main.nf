#!/usr/bin/env nextflow
nextflow.enable.dsl=2
projectId = params.projectId
read1AnalysisDataCode = params.read1AnalysisDataCode
read2AnalysisDataCode = params.read2AnalysisDataCode
referenceAnalysisDataCode = params.referenceAnalysisDataCode
pipelineId = params.pipelineId
pipelineCode = params.pipelineCode
userReference = params.userReference
storageSize = params.storageSize
fileUploadStatusCheckInterval = params.fileUploadStatusCheckInterval
analysisStatusCheckInterval = params.analysisStatusCheckInterval
readsFileUploadPath = params.readsFileUploadPath
readsPairFilesUploadPath = params.readsPairFilesUploadPath
referenceFileUploadPath = params.referenceFileUploadPath
localDownloadPath = params.localDownloadPath
icaUploadPath = params.icaUploadPath

process uploadFastqFilePairs {
    debug true
    input:
    tuple val(sampleId), file(reads)
    val(projectId)

    output:
    path "manifest.tsv", emit: manifestFile

    script:
    def (read_1_file, read_2_file) = reads
    """
    #!/bin/bash
    time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")
    read_1_file_id=""
    read_2_file_id=""
    ica_upload_path="/fastq/$sampleId/"

    printf "sampleId: $sampleId \n"
    printf "reads: $reads \n"
    printf "read_1_file: $read_1_file \n"
    printf "read_2_file: $read_2_file \n"

    read_1_file_response="read_1_file_response.txt"
    read_2_file_response="read_2_file_response.txt"

    touch \${read_1_file_response}
    touch \${read_2_file_response}

    printf "[\${time_stamp}]: "
    printf "Uploading read 1 file '${read_1_file}'... \n"
    read_1_upload_response=\$(icav2 projectdata upload ${read_1_file} \${ica_upload_path} --project-id ${projectId})
    echo "\${read_1_upload_response}" > \${read_1_file_response}

    printf "[\${time_stamp}]: "
    printf "Uploading read 2 file '${read_2_file}'... \n"
    read_2_upload_response=\$(icav2 projectdata upload ${read_2_file} \${ica_upload_path} --project-id ${projectId})
    echo "\${read_2_upload_response}" > \${read_2_file_response}

    # id of file starts with 'fil.'
    read_1_file_id=\$(cat \${read_1_file_response} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')
    read_2_file_id=\$(cat \${read_2_file_response} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')

    read_1_uploaded_file_data_response=\$(icav2 projectdata get \${read_1_file_id})
    if [[ \$? != 0 ]]; then
        printf "Failed to fetch data about file with id '\${read_1_file_id}'. \n"
        exit 1
    else
        read_1_uploaded_file_path=\$(echo \${read_1_uploaded_file_data_response} | jq -r ".details.path")
        printf "Path of uploaded file is '\${read_1_uploaded_file_path}'. \n"
    fi

    read_2_uploaded_file_data_response=\$(icav2 projectdata get \${read_2_file_id})
    if [[ \$? != 0 ]]; then
        printf "Failed to fetch data about file with id '\${read_2_file_id}'. \n"
        exit 1
    else
        read_2_uploaded_file_path=\$(echo \${read_2_uploaded_file_data_response} | jq -r ".details.path")
        printf "Path of uploaded file is '\${read_2_uploaded_file_path}'. \n"
    fi

    manifest_file="manifest.tsv"
    manifest_tab=\$(printf "\t")

    if ! [ -f \${manifest_file} ]; then
        echo "Manifest file does not exist. Creating one..."
        touch \${manifest_file}

        manifest_header_1="file_name"
        manifest_header_2="file_analysis_code"

        printf "[\${time_stamp}]: "
        printf "Writing file data to manifest...\n"

        printf "\${manifest_header_1} \${manifest_tab} \${manifest_header_2}\n" >> \${manifest_file}
        printf "\${read_1_file_id} \${manifest_tab} ${read1AnalysisDataCode}:\${read_1_file_id}\n" >> \${manifest_file}
        printf "\${read_1_file_id} \${manifest_tab} ${read2AnalysisDataCode}:\${read_2_file_id}\n" >> \${manifest_file}
    else
        printf "[\${time_stamp}]: "
        printf "Writing file data to existing manifest...\n"

        printf "\${read_1_file_id} \${manifest_tab} ${read1AnalysisDataCode}:\${read_1_file_id}\n" >> \${manifest_file}
        printf "\${read_1_file_id} \${manifest_tab} ${read2AnalysisDataCode}:\${read_2_file_id}\n" >> \${manifest_file}
    fi
    """
}

workflow {
    fastqFilePairs = Channel.fromFilePairs(readsPairFilesUploadPath, checkIfExists:true)
    fastqFilePairs.view()
    uploadFastqFilePairs(fastqFilePairs, params.projectId)

    referenceFilePath = Channel.fromPath(params.referenceFileUploadPath, checkIfExists: true)
}
