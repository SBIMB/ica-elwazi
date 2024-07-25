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
    maxForks 2
    input:
    tuple val(sampleId), file(reads)
    val(projectId)

    output:
    path "data.txt", emit: dataFile

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

    data_file="data.txt"

    if ! [ -f \${data_file} ]; then
        echo "Data file does not exist. Creating one..."
        touch \${data_file}
    fi

    printf "[\${time_stamp}]: "
    printf "Writing file data to existing data file...\n"

    printf "sampleId:${sampleId}\n" >> \${data_file}
    printf "${read1AnalysisDataCode}:\${read_1_file_id}\n" >> \${data_file}
    printf "${read2AnalysisDataCode}:\${read_2_file_id}\n" >> \${data_file}
    """
}

process uploadReferenceFile {
  debug true

  input:
  path(dataFile)
  path(referenceFileUploadPath)

  output:
  path "data.txt", emit: dataFile

  script:
  def reference_file = referenceFileUploadPath.baseName
  """
  #!/bin/bash
  time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")
  reference_file_id=${referenceFileId}
  reference_file_ica_path=${referenceFileIcaPath}

  get_reference_file_response_file="get_reference_file_response.txt"
  reference_file_upload_response_file="reference_file_upload_response.txt"

  touch \${reference_file_response}

  printf "[\${time_stamp}]: "
  printf "Checking if reference file id has been set in params.json...\n"
  if [ -z "${referenceFileId}" ]; then 
    printf "[\${time_stamp}]: "
    printf "Reference file id has been set in params.json. Getting reference file data JSON response...\n"
    get_reference_file_response=\$(icav2 projectdata get ${referenceFileId} --project-id ${projectId})
    echo "\${get_reference_file_response}" > \${get_reference_file_response_file}
  fi

  if grep -iq "No data found for path" \${get_reference_file_response_file}; then
    printf "[\${time_stamp}]: "
    printf "Reference file not found in ICA. Uploading reference file '${reference_file}'... \n"
    reference_file_upload_response=\$(icav2 projectdata upload ${referenceFileUploadPath} \${reference_file_ica_path} --project-id ${projectId})
    echo "\${reference_file_upload_response}" > \${reference_file_upload_response_file}

    printf "[\${time_stamp}]: "
    printf "Extracting file_id of reference file '${reference_file}' from upload response... \n"
    reference_file_id=\$(cat \${reference_file_upload_response_file} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')
  else
      printf "[\${time_stamp}]: "
      printf "Extracting file_id of reference file '${reference_file}' from get response... \n"
      reference_file_id=\$(cat \${get_reference_file_response_file} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')
  fi

  printf "[\${time_stamp}]: "
  printf "Writing file data to existing data file...\n"

  printf "${referenceAnalysisDataCode}:\${reference_file_id}\n" >> ${dataFile}
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
    path(dataFile)
    path(outputFolderId)

    output:
    stdout

    script:
    """
    sample_id=$(cat data_file.txt | grep -o 'sampleId:.*' | cut -f2- -d:)
    read_1_file_id=$(cat data_file.txt | grep -o 'read1:.*' | cut -f2- -d:)
    read_2_file_id=$(cat data_file.txt | grep -o 'read2:.*' | cut -f2- -d:)
    reference_file_id=$(cat data_file.txt | grep -o 'ref_tar:.*' | cut -f2- -d:)

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting uploaded read 1 file with ID '\${read_1_file_id}'...\n"
    icav2 projectdata delete \${read_1_file_id}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting uploaded read 2 file with ID '\${read_2_file_id}'...\n"
    icav2 projectdata delete \${read_2_file_id}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting uploaded reference file with ID '\${reference_file_id}'...\n"
    icav2 projectdata delete \${reference_file_id}

    folderId=\$(cat ${outputFolderId})
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting analysis output folder with ID '\${folderId}'...\n"
    icav2 projectdata delete \${folderId}

    printf "Uploaded file and analysis output folder successfully deleted.\n"
    """
}

workflow {
    fastqFilePairs = Channel.fromFilePairs(readsPairFilesUploadPath, checkIfExists:true)
    referenceFilePath = Channel.fromPath(params.referenceFileUploadPath, checkIfExists: true)

    uploadFastqFilePairs(fastqFilePairs, params.projectId)
    uploadReferenceFile(uploadFastqFilePairs.out.dataFile, referenceFilePath)
    startAnalysis(uploadReferenceFile.out.dataFile)
    checkAnalysisStatus(startAnalysis.out.analysisResponse, params.analysisStatusCheckInterval)
    downloadAnalysisOutput(checkAnalysisStatus.out.analysisOutputFolderId, params.localDownloadPath)
    deleteData(uploadFile.out.fileUploadResponse.view(), downloadAnalysisOutput.out.outputFolderId)
}
