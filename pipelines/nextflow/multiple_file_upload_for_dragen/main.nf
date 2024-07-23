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
referenceFileUploadPath = params.referenceFileUploadPath
localDownloadPath = params.localDownloadPath
icaUploadPath = params.icaUploadPath

process uploadFiles {
    debug true
    input:
    path(readsFilePath)
    path(referenceFilePath)
    val(projectId)

    output:
    path "manifest.tsv", emit: manifestFile

    script:
    referenceFileName = referenceFilePath.baseName
    """
    #!/bin/bash
    time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")
    read_1_file=""
    read_2_file=""

    for file in "${readsFilePath}"/*; do
        if [[ "\${file}" == *"_R1_"* ]]; then
            read_1_file=\${file}
            printf "Read 1 file is '\${read_1_file}'.\n"
        elif [[ "\${file}" == *"_R2_"* ]]; then
            read_2_file=\${file}
            printf "Read 2 file is '\${read_2_file}'.\n"
        else
            printf "No read files present in directory.\n"
        fi
    done

    read_1_file_response="read_1_file_response.txt"
    read_2_file_response="read_2_file_response.txt"
    reference_file_response="reference_file_response.txt"

    touch \${read_1_file_response}
    touch \${read_2_file_response}
    touch \${reference_file_response}

    printf "[\${time_stamp}]: "
    printf "Uploading read 1 file '\${read_1_file}'... \n"
    read_1_upload_response=\$(icav2 projectdata upload \${read_1_file} ${icaUploadPath} --project-id ${projectId})
    echo "\${read_1_upload_response}" > \${read_1_file_response}

    printf "[\${time_stamp}]: "
    printf "Uploading read 2 file '\${read_2_file}'... \n"
    read_2_upload_response=\$(icav2 projectdata upload \${read_2_file} ${icaUploadPath} --project-id ${projectId})
    echo "\${read_2_upload_response}" > \${read_2_file_response}

    printf "[\${time_stamp}]: "
    printf "Uploading reference file '${referenceFileName}'... \n"
    reference_file_upload_response=\$(icav2 projectdata upload ${referenceFilePath} ${icaUploadPath} --project-id ${projectId})
    echo "\${reference_file_upload_response}" > \${reference_file_response}.txt

    # id of file starts with 'fil.'
    read_1_file_id=\$(cat \${read_1_file_response} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')
    read_2_file_id=\$(cat \${read_2_file_response} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')
    reference_file_id=\$(cat \${read_1_file_response} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')

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

    uploaded_reference_file_data_response=\$(icav2 projectdata get \${reference_file_id})
    if [[ \$? != 0 ]]; then
        printf "Failed to fetch data about reference file with id '\${reference_file_id}'. \n"
        exit 1
    else
        reference_file_upload_path=\$(echo \${uploaded_reference_file_data_response} | jq -r ".details.path")
        printf "Path of uploaded reference file is '\${reference_file_upload_path}'. \n"
    fi

    manifest_file="manifest.tsv"
    manifest_tab=\$(printf "\t")

    if ! [ -f \${manifest_file} ]; then
        echo "Manifest file does not exist. Creating one..."
        touch \${manifest_file}

        manifest_header_1="file_name"
        manifest_header_2="file_id"
        manifest_header_3="file_ica_path"
        manifest_header_4="file_analysis_code"

        printf "[\${time_stamp}]: "
        printf "Writing file data to manifest...\n"

        printf "\${manifest_header_1} \${manifest_tab} \${manifest_header_2} \${manifest_tab} \${manifest_header_3} \${manifest_tab} \${manifest_header_4}\n" >> \${manifest_file}
        printf "\${read_1_file} \${manifest_tab} \${read_1_file_id} \${manifest_tab} \${read_1_uploaded_file_path} \${manifest_tab} ${read1AnalysisDataCode}:\${read_1_file_id}\n" >> \${manifest_file}
        printf "\${read_2_file} \${manifest_tab} \${read_2_file_id} \${manifest_tab} \${read_2_uploaded_file_path} \${manifest_tab} ${read2AnalysisDataCode}:\${read_2_file_id}\n" >> \${manifest_file}
        printf "${referenceFileName} \${manifest_tab} \${reference_file_id} \${manifest_tab} \${reference_file_upload_path} \${manifest_tab} ${referenceAnalysisDataCode}:\${reference_file_id}\n" >> \${manifest_file}
    else
        printf "[\${time_stamp}]: "
        printf "Writing file data to existing manifest...\n"

        printf "\${read_1_file} \${manifest_tab} \${read_1_file_id} \${manifest_tab} \${read_1_uploaded_file_path} \${manifest_tab} ${read1AnalysisDataCode}:\${read_1_file_id}\n" >> \${manifest_file}
        printf "\${read_2_file} \${manifest_tab} \${read_2_file_id} \${manifest_tab} \${read_2_uploaded_file_path} \${manifest_tab} ${read2AnalysisDataCode}:\${read_2_file_id}\n" >> \${manifest_file}
        printf "${referenceFileName} \${manifest_tab} \${reference_file_id} \${manifest_tab} \${reference_file_upload_path} \${manifest_tab} ${referenceAnalysisDataCode}:\${reference_file_id}\n" >> \${manifest_file}
    fi

    """
}

process constructFileReference {
    debug true
    
    input:
    path(fileUploadResponse)
    val(analysisDataCode)

    output:
    path "fileReference.txt", emit: fileRef

    script:
    fileReference = ""
    """
    #!/bin/bash
    
    fileUploadStatusCheckCount=0
    fileUploadStatusCheckLimit=10
    fileUploadStatus="PARTIAL"

    fileId=\$(cat ${fileUploadResponse} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\${timeStamp}]: Checking status of uploaded file with id '\${fileId}' every ${fileUploadStatusCheckInterval} seconds, until status is 'AVAILABLE'..."
    while true;
    do
        ((\${fileUploadStatusCheckCount}+=1))
        uploadedFileResponse=\$(icav2 projectdata get \${fileId})

        echo "Checking status of file with id '\${fileId}'..."
        fileUploadStatus=\$(echo \${uploadedFileResponse} | jq -r ".details.status")

        echo "Current status of uploaded file is '\${fileUploadStatus}'."
        if [[ \${fileUploadStatus} == "AVAILABLE" ]]; then
            echo "Uploaded file is AVAILABLE"
            fileReference="${analysisDataCode}:\${fileId}"
            echo "File Reference:\${fileReference}"
            touch fileReference.txt
            echo "\${fileReference}" > fileReference.txt
            break;
        elif [[ \${fileUploadStatusCheckCount} -gt \${fileUploadStatusCheckLimit} ]]; then
            echo "Uploaded file status has been checked more than \${fileUploadStatusCheckLimit} times. Stopping..."
            break;
        else
            echo "Uploaded file is still not AVAILABLE. Checking again..."
        fi
        sleep ${fileUploadStatusCheckInterval};
    done
    """
}

process startAnalysis {
    debug true
    
    input:
    path(fileRef)

    output:
    path "analysisResponse.txt", emit: analysisResponse

    script:
    analysisResponse = ""

    """
    #!/bin/bash
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\${timeStamp}]: Starting Nextflow analysis..."

    fileReference=\$(cat ${fileRef})

    echo "File Ref: \${fileReference}"

    analysisResponse=\$(icav2 projectpipelines start nextflow ${pipelineId} \
        --user-reference ${userReference} \
        --project-id ${projectId} \
        --storage-size ${storageSize} \
        --input \${fileReference})

    echo "\${analysisResponse}" > analysisResponse.txt
    """
}

process checkAnalysisStatus {
    debug true
    
    input:
    path(analysisResponse)
    val(analysisStatusCheckInterval)

    output:
    path "analysisOutputFolderId.txt", emit: analysisOutputFolderId

    script:
    analysisOutputFolderId = ""
    """
    #!/bin/bash

    analysisStatusCheckCount=0
    analysisStatusCheckLimit=10
    analysisStatus="REQUESTED"

    analysisId=\$(cat ${analysisResponse} | jq -r ".id")
    analysisRef=\$(cat ${analysisResponse} | jq -r ".reference")
    
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\${timeStamp}]: Checking status of analysis with id '\${analysisId}' every ${analysisStatusCheckInterval} seconds, until status is 'SUCCEEDED'..."
    while true;
    do
        ((\${StatusCheckCount}+=1))
        updatedAnalysisResponse=\$(icav2 projectanalyses get \${analysisId})

        echo "Checking status of analysis with reference '\${analysisRef}'..."
        analysisStatus=\$(echo \${updatedAnalysisResponse} | jq -r ".status")

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        echo "[\${timeStamp}]: Current status of analysis is '\${analysisStatus}'..."

        if [[ \${analysisStatus} == "SUCCEEDED" ]]; then
            echo "Analysis SUCCEEDED"
            echo "Fetching analysis output response..."
            analysisOutputResponse=\$(icav2 projectanalyses output \$analysisId)
            analysisOutputFolderId=\$(echo \${analysisOutputResponse} | jq -r ".items[].data[].dataId")
            echo "Analysis output folder ID is '\${analysisOutputFolderId}'"

            touch analysisOutputFolderId.txt
            echo "\${analysisOutputFolderId}" > analysisOutputFolderId.txt
            break;

        elif [[ \${analysisStatus} == "FAILED" ]]; then
            echo "Analysis FAILED \n"
            break;

        elif [[ \${analysisStatus} == "FAILED_FINAL" ]]; then
            echo "Analysis FAILED_FINAL"
            break;

        elif [[ \${analysisStatus} == "ABORTED" ]]; then
            echo "Analysis ABORTED"
            break;

        elif [[ \${analysisStatusCheckCount} -gt \${analysisStatusCheckLimit} ]]; then
            echo "Analysis status has been checked more than \${analysisStatusCheckLimit} times. Stopping..."
            break;

        else
            echo "Analysis still in progress..."
        fi

        sleep ${analysisStatusCheckInterval};
    done
    """
}

process downloadAnalysisOutput {
    debug true
    
    input:
    path(analysisOutputFolderId)
    val(localDownloadPath)

    output:
    path "outputFolderId.txt", emit: outputFolderId

    script:
    outputFolderId = ""
    """
    #!/bin/bash

    outputFolderId=\$(cat ${analysisOutputFolderId})

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\${timeStamp}]: Downloading analysis output folder with ID '\${outputFolderId}' to '${localDownloadPath}'..."

    icav2 projectdata download \${outputFolderId} ${localDownloadPath}

    touch outputFolderId.txt
    echo "\${outputFolderId}" > outputFolderId.txt
    """
}

process deleteData {
    debug true

    input:
    path(fileUploadResponse)
    path(outputFolderId)

    output:
    stdout

    script:
    """
    fileId=\$(cat ${fileUploadResponse} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\${timeStamp}]: Deleting uploaded file with ID '\${fileId}'..."
    icav2 projectdata delete \${fileId}

    folderId=\$(cat ${outputFolderId})
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\${timeStamp}]: Deleting analysis output folder with ID '\${folderId}'..."
    icav2 projectdata delete \${folderId}

    echo "Uploaded file and analysis output folder successfully deleted."
    """
}

workflow {
    readsFilePath = Channel.fromPath(params.readsFileUploadPath, checkIfExists: true)
    referenceFilePath = Channel.fromPath(params.referenceFileUploadPath, checkIfExists: true)

    uploadFiles(readsFilePath, referenceFilePath, params.projectId)
    
    // constructFileReference(uploadFiles.out.fileUploadResponse.view(), params.analysisDataCode)

    // startAnalysis(constructFileReference.out.fileRef)

    // checkAnalysisStatus(startAnalysis.out.analysisResponse, params.analysisStatusCheckInterval)

    // downloadAnalysisOutput(checkAnalysisStatus.out.analysisOutputFolderId, params.localDownloadPath)

    // deleteData(uploadFiles.out.fileUploadResponse.view(), downloadAnalysisOutput.out.outputFolderId)
}
