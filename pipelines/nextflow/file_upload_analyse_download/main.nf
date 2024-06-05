#!/usr/bin/env nextflow
nextflow.enable.dsl=2
filePath = Channel.fromPath("NZ_GG704948.fa", checkIfExists: true)
projectId = params.projectId
analysisDataCode = params.analysisDataCode
pipelineId = params.pipelineId
pipelineCode = params.pipelineCode
userReference = params.userReference
storageSize = params.storageSize
fileUploadStatusCheckInterval = params.fileUploadStatusCheckInterval
analysisStatusCheckInterval = params.analysisStatusCheckInterval
localDownloadPath = params.localDownloadPath

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
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    echo "[\${timeStamp}]: Uploading file '${fileName}' to project with id '${projectId}'..."

    touch ${fileName}.txt

    fileUploadResponse=\$(icav2 projectdata upload ${filePath} --project-id ${projectId})

    echo "Successfully uploaded file '${fileName}' to project with id '${projectId}'."

    echo "\${fileUploadResponse}" > ${fileName}.txt
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
        echo "Current status of analysis is '\${analysisStatus}'..."

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
    uploadFile(filePath, params.projectId)
    
    constructFileReference(uploadFile.out.fileUploadResponse.view(), params.analysisDataCode)

    startAnalysis(constructFileReference.out.fileRef)

    checkAnalysisStatus(startAnalysis.out.analysisResponse, params.analysisStatusCheckInterval)

    downloadAnalysisOutput(checkAnalysisStatus.out.analysisOutputFolderId, params.localDownloadPath)

    deleteData(uploadFile.out.fileUploadResponse.view(), downloadAnalysisOutput.out.outputFolderId)
}
