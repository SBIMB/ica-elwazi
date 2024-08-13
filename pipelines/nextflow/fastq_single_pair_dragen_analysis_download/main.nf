#!/usr/bin/env nextflow
nextflow.enable.dsl=2
projectId = params.projectId
pipelineId = params.pipelineId
pipelineCode = params.pipelineCode
userReference = params.userReference
storageSize = params.storageSize
fileUploadStatusCheckInterval = params.fileUploadStatusCheckInterval
analysisStatusCheckInterval = params.analysisStatusCheckInterval
analysisStatusCheckLimit = params.analysisStatusCheckLimit
read1FileId = params.read1FileId
read2FileId = params.read2FileId
analysisId = params.analysisId
readsFileUploadPath = params.readsFileUploadPath
referenceFileId = params.referenceFileId
readsPairFilesUploadPath = params.readsPairFilesUploadPath
referenceFileUploadPath = params.referenceFileUploadPath
localDownloadPath = params.localDownloadPath

process createDataFile {
    debug true

    output:
    path "data.txt", emit: dataFile

    script:
    """
    #!/bin/bash
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Creating data file...\n"
    
    data_file="data.txt"

    if ! [ -f \${data_file} ]; then
        echo "Data file does not exist. Creating one..."
        touch \${data_file}
    fi

    printf "[\${time_stamp}]: "
    printf "Writing file data to existing data file...\n"

    printf "read1:${read1FileId}\n" >> \${data_file}
    printf "read2:${read2FileId}\n" >> \${data_file}
    printf "analysisId:${analysisId}\n" >> \${data_file}
    """

}
process checkAnalysisStatus {
    debug true
    
    input:
    path(dataFile)

    output:
    path "data.txt", emit: dataFile

    script:
    """
    #!/bin/bash

    analysis_status_check_count=0
    analysis_status="REQUESTED"

    analysis_id=\$(cat ${dataFile} | grep -o 'analysisId:.*' | cut -f2- -d:)
    echo "\${analysis_id}"
    while_loop_completion_message="Exiting WHILE loop. Moving on to 'downloadAnalysisOutput' process..."

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Checking status of analysis with id '\${analysis_id}' every ${analysisStatusCheckInterval} seconds, until status is 'SUCCEEDED'...\n"

    while true;
    do
        ((analysis_status_check_count +=1 ))
        updated_analysis_response=\$(icav2 projectanalyses get \${analysis_id})
        analysis_status=\$(echo \${updated_analysis_response} | jq -r ".status")

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Current status of analysis is '\${analysis_status}'...\n"

        if [[ \${analysis_status} == "SUCCEEDED" ]]; then
            printf "Analysis SUCCEEDED\n"
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

        elif [[ \${analysis_status_check_count} -gt ${analysisStatusCheckLimit} ]]; then
            printf "Analysis status has been checked more than ${analysisStatusCheckLimit} times. Stopping...\n"
            printf "analysisStatus:TIMEOUT\n" >> ${dataFile}
            break;

        else
            printf "Analysis still in progress...\n"
        fi

        sleep ${analysisStatusCheckInterval};
    done

    printf "analysisStatus:\${analysis_status}\n" >> ${dataFile}

    printf "\${while_loop_completion_message}\n"
    """
}

process downloadAnalysisOutput {
    debug true
    
    input:
    path(dataFile)

    output:
    path "data.txt", emit: dataFile

    script:
    """
    #!/bin/bash
    analysis_id=\$(cat ${dataFile} | grep -o 'analysisId:.*' | cut -f2- -d:)

    printf "[\${time_stamp}]: "
    printf "Fetching analysis output response...\n"
    analysis_output_response=\$(icav2 projectanalyses output \${analysis_id})
    analysis_output_folder_id=\$(echo \${analysis_output_response} | jq -r ".items[].data[].dataId")
    printf "Analysis output folder ID is '\${analysis_output_folder_id}'\n"
    printf "Writing id of analysis output folder to existing data file...\n"
    printf "outputFolderId:\${analysis_output_folder_id}\n" >> ${dataFile}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Downloading analysis output folder with ID '\${analysis_output_folder_id}' to '${localDownloadPath}'...\n"

    icav2 projectdata download \${analysis_output_folder_id} ${localDownloadPath}

    printf "downloadComplete:true\n" >> ${dataFile}
    """
}

process deleteData {
    debug true

    input:
    path(dataFile)

    output:
    stdout

    script:
    """
    read_1_file_id=\$(cat ${dataFile} | grep -o 'read1:.*' | cut -f2- -d:)
    read_2_file_id=\$(cat ${dataFile} | grep -o 'read2:.*' | cut -f2- -d:)
    analysis_output_folder_id=\$(cat ${dataFile} | grep -o 'outputFolderId:.*' | cut -f2- -d:)
    download_complete=\$(cat ${dataFile} | grep -o 'downloadComplete:.*' | cut -f2- -d:)

    if [ "\${download_complete}" = "true" ]; then  
        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Deleting uploaded read 1 file with ID '\${read_1_file_id}'...\n"
        icav2 projectdata delete \${read_1_file_id}
        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Deleting uploaded read 2 file with ID '\${read_2_file_id}'...\n"
        icav2 projectdata delete \${read_2_file_id}

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Deleting analysis output folder with ID '\${analysis_output_folder_id}'...\n"
        icav2 projectdata delete \${analysis_output_folder_id}
    fi
    printf "Uploaded file and analysis output folder successfully deleted.\n"
    """
}

workflow {
    // dataFilePath = Channel.fromPath("${projectDir}/data.txt", checkIfExists: true)
    createDataFile()
    // def dataFileChannel = createDataFile.out.dataFile
    def dataFileChannel = Channel.fromPath("${projectDir}/data.txt")
    dataFileChannel.view()

    println "dataFileChannel: ${dataFileChannel}"

    checkAnalysisStatus(createDataFile.out.dataFile)

    downloadAnalysisOutput(checkAnalysisStatus.out.dataFile)
    deleteData(downloadAnalysisOutput.out.dataFile)
}
