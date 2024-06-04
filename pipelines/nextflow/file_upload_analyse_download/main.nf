#!/usr/bin/env nextflow
nextflow.enable.dsl=2
filePath = Channel.fromPath("NZ_GG704945.fa", checkIfExists: true)
projectId = params.projectId
analysisDataCode = params.analysisDataCode
pipelineId = params.pipelineId
pipelineCode = params.pipelineCode
userReference = params.userReference
storageSize = params.storageSize
analysisStatusCheckInterval = params.analysisStatusCheckInterval

process uploadFile {
    debug true
    input:
    path(filePath)
    val(projectId)

    output:
    path "${fileName}.txt", emit: fileUploadResponse

    script:
    fileName = filePath.baseName
    """
    #!/bin/bash

    echo "Uploading file '${fileName}' to project with id '${projectId}'..."

    touch ${fileName}.txt

    fileUploadResponse=\$(icav2 projectdata upload ${filePath} --project-id ${projectId})

    echo "Successfully uploaded file '${fileName}' to project with id '${projectId}'."

    echo "\${fileUploadResponse}" > ${fileName}.txt
    """
}

process checkFileUploadStatus {
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
    
    fileId=\$(cat ${fileUploadResponse} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')

    fileReference="${analysisDataCode}:\${fileId}"

    echo "File Reference:\${fileReference}"

    touch fileReference.txt

    echo "\${fileReference}" > fileReference.txt
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
    echo "Starting Nextflow analysis..."

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
    val(analysisResponse)
    val(analysisStatusCheckInterval)

    output:
    val analysisOutputFolderId, emit: analysisOutputFolderId

    script:
    analysisOutputFolderId = ""
    """
    #!/bin/bash

    analysisStatusCheckCount=0
    analysisStatusCheckLimit=10
    analysisStatus="REQUESTED"

    analysisId=\$(cat ${analysisResponse} | jq -r ".id")
    analysisRef=\$(cat ${analysisResponse} | jq -r ".reference")
    
    echo "Checking status of analysis with id '\${analysisId}' every ${analysisStatusCheckInterval} seconds, until status is 'SUCCEEDED'..."
    while true;
    do
        ((\${StatusCheckCount}+=1))
        updatedAnalysisResponse=\$(icav2 projectanalyses get \${analysisId})

        echo "Checking status of analysis with reference '\${analysisRef}'..."
        analysisStatus=\$(echo \${updatedAnalysisResponse} | jq -r ".items[] | select(.reference == "\${analysisRef}").status")
        echo "Current status of analysis is '\${analysisStatus}'..."

        if [[ \${analysisStatus} == "SUCCEEDED" ]]; then
            echo "Analysis SUCCEEDED"
            echo "Fetching analysis output response..."
            analysisOutputResponse=\$(icav2 projectanalyses output \$analysisId)
            analysisOutputFolderId=\$(echo \${analysisOutputResponse} | jq -r ".items[].data[].dataId")
            echo "Analysis output folder ID is '\${analysisOutputFolderId}'"
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
    val(analysisOutputFolderId)
    val(localDownloadPath)

    output:
    stdout

    script:
    """
    #!/bin/bash

    echo "Downloading analysis output folder with ID '${analysisOutputFolderId}' to '${localDownloadPath}'..."

    icav2 projectdata download ${analysisOutputFolderId} ${local_download_path}
    """
}

workflow {
    uploadFile(filePath, params.projectId)
    
    checkFileUploadStatus(uploadFile.out.fileUploadResponse.view(), params.analysisDataCode)

    startAnalysis(checkFileUploadStatus.out.fileRef)

    checkAnalysisStatus(startAnalysis.out.analysisResponse.view(), params.analysisStatusCheckInterval)

    downloadAnalysisOutput(checkAnalysisStatus.out.analysisOutputFolderId.view(), params.localDownloadPath)
}
