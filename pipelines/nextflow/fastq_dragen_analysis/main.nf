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
sampleId = params.sampleId
read1FileId = params.read1FileId
read2FileId = params.read2FileId
referenceFileId = params.referenceFileId
localDownloadPath = params.localDownloadPath
icaUploadPath = params.icaUploadPath

process checkFileStatus {
    debug true
    
    output:
    path "data.txt", emit: dataFile

    script:

    """
    #!/bin/bash
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Getting status of uploaded files...\n"

    read_1_file_data_response=\$(icav2 projectdata get ${read1FileId})
    read_1_file_status=\$(echo \${read_1_file_data_response} | jq -r ".details.status")

    read_2_file_data_response=\$(icav2 projectdata get ${read2FileId})
    read_2_file_status=\$(echo \${read_2_file_data_response} | jq -r ".details.status")

    reference_file_data_response=\$(icav2 projectdata get ${referenceFileId})
    reference_file_status=\$(echo \${reference_file_data_response} | jq -r ".details.status")

    if [[ \${read_1_file_status} == "AVAILABLE" ]]; then
        printf "Read 1 file is AVAILABLE\n"
        read_1_analysis_code="${read1AnalysisDataCode}:${read1FileId}"
    else
        printf "Read 1 file is not AVAILABLE\n"
    fi

    if [[ \${read_2_file_status} == "AVAILABLE" ]]; then
        printf "Read 2 file is AVAILABLE\n"
        read_2_analysis_code="${read2AnalysisDataCode}:${read2FileId}"
    else
        printf "Read 2 file is not AVAILABLE\n"
    fi

    if [[ \${reference_file_status} == "AVAILABLE" ]]; then
        printf "Reference file is AVAILABLE\n"
        reference_file_analysis_code="${referenceAnalysisDataCode}:${referenceFileId}"
    else
        printf "Reference file is not AVAILABLE\n"
    fi
    
    data_file="data.txt"
    touch \${data_file}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Writing data to data file...\n"

    printf "sampleId:${sampleId}\n" >> \${data_file}
    printf "\${read_1_analysis_code}\n" >> \${data_file}
    printf "\${read_2_analysis_code}\n" >> \${data_file}
    printf "\${reference_file_analysis_code}\n" >> \${data_file}
    """
}

process startAnalysis {
    debug true
    
    input:
    path(dataFile)

    output:
    path "analysisResponse.txt", emit: analysisResponse

    script:
    analysisResponse = ""

    """
    #!/bin/bash

    sample_id=\$(cat ${dataFile} | grep -E "sampleId")

    read1_analysis_code=\$(cat ${dataFile} | grep -E "read1")
    read2_analysis_code=\$(cat ${dataFile} | grep -E "read2")
    reference_analysis_code=\$(cat ${dataFile} | grep -E "ref_tar")

    output_directory="/output/${sampleId}/"

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Starting Nextflow analysis...\n"
    analysisResponse=\$(icav2 projectpipelines start nextflow ${pipelineId} \
        --user-reference ${userReference} \
        --project-id ${projectId} \
        --storage-size ${storageSize} \
        --input \${read1_analysis_code} \
        --input \${read2_analysis_code} \
        --input \${reference_analysis_code} \
        --parameters enable-variant-caller:true \
        --parameters RGID:Illumina_RGID \
        --parameters RGSM:${sampleId} \
        --parameters output-directory:\${output_directory} \
        --parameters output-file-prefix:${sampleId}) 

    touch "analysisResponse.txt"
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

    analysis_status_check_count=0
    analysis_status_check_limit=10
    analysis_status="REQUESTED"

    analysis_id=\$(cat ${analysisResponse} | jq -r ".id")
    analysis_ref=\$(cat ${analysisResponse} | jq -r ".reference")
    
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Checking status of analysis with id '\${analysis_id}' every ${analysisStatusCheckInterval} seconds, until status is 'SUCCEEDED'...\n"
    while true;
    do
        ((\${StatusCheckCount}+=1))
        updatedAnalysisResponse=\$(icav2 projectanalyses get \${analysis_id})

        printf "Checking status of analysis with reference '\${analysis_ref}'...\n"
        analysis_status=\$(echo \${updatedAnalysisResponse} | jq -r ".status")

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Current status of analysis is '\${analysis_status}'...\n"

        if [[ \${analysis_status} == "SUCCEEDED" ]]; then
            printf "Analysis SUCCEEDED\n"
            printf "Fetching analysis output response...\n"
            analysisOutputResponse=\$(icav2 projectanalyses output \$analysis_id)
            analysisOutputFolderId=\$(echo \${analysisOutputResponse} | jq -r ".items[].data[].dataId")
            printf "Analysis output folder ID is '\${analysisOutputFolderId}'\n"

            touch analysisOutputFolderId.txt
            echo "\${analysisOutputFolderId}" > analysisOutputFolderId.txt
            break;

        elif [[ \${analysis_status} == "FAILED" ]]; then
            printf "Analysis FAILED \n"
            break;

        elif [[ \${analysis_status} == "FAILED_FINAL" ]]; then
            printf "Analysis FAILED_FINAL\n"
            break;

        elif [[ \${analysis_status} == "ABORTED" ]]; then
            printf "Analysis ABORTED\n"
            break;

        elif [[ \${analysis_status_check_count} -gt \${analysis_status_check_limit} ]]; then
            printf "Analysis status has been checked more than \${analysis_status_check_limit} times. Stopping...\n"
            break;

        else
            printf "Analysis still in progress...\n"
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
    printf "[\${timeStamp}]: Downloading analysis output folder with ID '\${outputFolderId}' to '${localDownloadPath}'...\n"

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
    printf "[\${timeStamp}]: Deleting uploaded file with ID '\${fileId}'...\n"
    icav2 projectdata delete \${fileId}

    folderId=\$(cat ${outputFolderId})
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting analysis output folder with ID '\${folderId}'...\n"
    icav2 projectdata delete \${folderId}

    printf "Uploaded file and analysis output folder successfully deleted.\n"
    """
}


workflow {
    checkFileStatus()
    startAnalysis(checkFileStatus.out.dataFile)
    checkAnalysisStatus(startAnalysis.out.analysisResponse, params.analysisStatusCheckInterval)
    downloadAnalysisOutput(checkAnalysisStatus.out.analysisOutputFolderId, params.localDownloadPath)
    // deleteData(uploadFile.out.fileUploadResponse.view(), downloadAnalysisOutput.out.outputFolderId)
}
