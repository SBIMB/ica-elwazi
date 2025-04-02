#!/usr/bin/env nextflow
nextflow.enable.dsl=2
projectId = params.projectId
read1AnalysisDataCode = params.read1AnalysisDataCode
read2AnalysisDataCode = params.read2AnalysisDataCode
fastqsAnalysisDataCode = params.fastqsAnalysisDataCode
fastqListDataCode = params.fastqListDataCode
referenceAnalysisDataCode = params.referenceAnalysisDataCode
pipelineId = params.pipelineId
pipelineCode = params.pipelineCode
userReference = params.userReference
storageSize = params.storageSize
hashTableConfigFile = params.hashTableConfigFile
referenceDirectory = params.referenceDirectory
intermediateResultsDirectory = params.intermediateResultsDirectory
fileUploadStatusCheckInterval = params.fileUploadStatusCheckInterval
analysisStatusCheckInterval = params.analysisStatusCheckInterval
analysisStatusCheckLimit = params.analysisStatusCheckLimit
readsFileUploadPath = params.readsFileUploadPath
referenceFileId = params.referenceFileId
readsPairFilesUploadPath = params.readsPairFilesUploadPath
localDownloadPath = params.localDownloadPath

process uploadFastqFilePairs {
    debug true
    maxForks 4
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
        read_1_uploaded_file_name=\$(echo \${read_1_uploaded_file_data_response} | jq -r ".details.name")
        printf "Name of uploaded file is '\${read_1_uploaded_file_name}'. \n"
        read_1_uploaded_file_path=\$(echo \${read_1_uploaded_file_data_response} | jq -r ".details.path")
        printf "Path of uploaded file is '\${read_1_uploaded_file_path}'. \n"
    fi

    read_2_uploaded_file_data_response=\$(icav2 projectdata get \${read_2_file_id})
    if [[ \$? != 0 ]]; then
        printf "Failed to fetch data about file with id '\${read_2_file_id}'. \n"
        exit 1
    else
        read_2_uploaded_file_name=\$(echo \${read_2_uploaded_file_data_response} | jq -r ".details.name")
        printf "Name of uploaded file is '\${read_2_uploaded_file_name}'. \n"
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
    printf "read1Name:\${read_1_uploaded_file_name}\n" >> \${data_file}
    printf "read2Name:\${read_2_uploaded_file_name}\n" >> \${data_file}
    """
}

process uploadFastqFileList {
    debug true
    input:
    path(dataFile)

    output:
    path "data.txt", emit: dataFile

    script:
    """
    #!/bin/bash
    time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")
    fastq_list_file_id=""

    printf "[\${time_stamp}]: "
    printf "Creating FASTQ list file... \n"
    sample_id=\$(cat ${dataFile} | grep -o 'sampleId:.*' | cut -f2- -d:)
    read_1_file_name=\$(cat ${dataFile} | grep -o 'read1Name:.*' | cut -f2- -d:)
    read_2_file_name=\$(cat ${dataFile} | grep -o 'read2Name:.*' | cut -f2- -d:)

    fastq_list_file="fastq-list-\${sample_id}.csv"
    touch \${fastq_list_file}

    printf "[\${time_stamp}]: "
    printf "Writing to FASTQ list file... \n"
    printf "RGID,RGSM,RGLB,Lane,Read1File,Read2File\n" >> \${fastq_list_file}
    printf "\${sample_id},\${sample_id},RGLB,1,\${read_1_file_name},\${read_2_file_name}\n" >> \${fastq_list_file}

    fastq_list_file_response="fastq_list_file_response.txt"
    touch \${fastq_list_file_response}

    printf "[\${time_stamp}]: "
    printf "Uploading FASTQ list file... \n"
    fastq_list_upload_response=\$(icav2 projectdata upload \${fastq_list_file} --project-id ${projectId})
    echo "\${fastq_list_upload_response}" > \${fastq_list_file_response}

    # id of file starts with 'fil.'
    fastq_list_file_id=\$(cat \${fastq_list_file_response} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')

    fastq_list_uploaded_file_data_response=\$(icav2 projectdata get \${fastq_list_file_id})
    if [[ \$? != 0 ]]; then
        printf "Failed to fetch data about FASTQ list file with id '\${fastq_list_file_id}'. \n"
        exit 1
    fi

    printf "[\${time_stamp}]: "
    printf "Writing FASTQ list file id to existing data file...\n"

    printf "${fastqListDataCode}:\${fastq_list_file_id}\n" >> ${dataFile}
    """
}

process getReferenceFile {
  debug true

  input:
  path(dataFile)

  output:
  path "data.txt", emit: dataFile

  script:
  def reference_file_name = ""
  """
  #!/bin/bash
  time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")

  get_reference_file_response_file="get_reference_file_response.txt"
  touch \${get_reference_file_response_file}

  printf "[\${time_stamp}]: "
  printf "Getting reference file data JSON response...\n"
  get_reference_file_response=\$(icav2 projectdata get ${referenceFileId} --project-id ${projectId})
  reference_file_name=\$(echo \${get_reference_file_response} | jq -r ".details.name")
  echo "\${get_reference_file_response}" > \${get_reference_file_response_file}

  printf "[\${time_stamp}]: "
  printf "Writing reference file '\${reference_file_name}' id to existing data file...\n"

  printf "${referenceAnalysisDataCode}:${referenceFileId}\n" >> ${dataFile}
  """
}

process checkFileStatus {
    debug true
    
    input:
    path(dataFile)

    output:
    path "data.txt", emit: dataFile

    script:
    """
    #!/bin/bash
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Getting status of uploaded files...\n"

    read_1_file_id=\$(cat ${dataFile} | grep -o 'read1:.*' | cut -f2- -d:)
    read_2_file_id=\$(cat ${dataFile} | grep -o 'read2:.*' | cut -f2- -d:)
    fastq_list_file_id=\$(cat ${dataFile} | grep -o 'fastq_list:.*' | cut -f2- -d:)
    reference_file_id=\$(cat ${dataFile} | grep -o 'ref_tar:.*' | cut -f2- -d:)

    read_1_file_data_response=\$(icav2 projectdata get \${read_1_file_id})
    read_1_file_status=\$(echo \${read_1_file_data_response} | jq -r ".details.status")

    read_2_file_data_response=\$(icav2 projectdata get \${read_2_file_id})
    read_2_file_status=\$(echo \${read_2_file_data_response} | jq -r ".details.status")

    fastq_list_file_data_response=\$(icav2 projectdata get \${fastq_list_file_id})
    fastq_list_file_status=\$(echo \${fastq_list_file_data_response} | jq -r ".details.status")

    reference_file_data_response=\$(icav2 projectdata get \${reference_file_id})
    reference_file_status=\$(echo \${reference_file_data_response} | jq -r ".details.status")

    if [[ \${read_1_file_status} == "AVAILABLE" ]]; then
        printf "Read 1 file is AVAILABLE\n"
    else
        printf "Read 1 file is not AVAILABLE\n"
        exit 1
    fi

    if [[ \${read_2_file_status} == "AVAILABLE" ]]; then
        printf "Read 2 file is AVAILABLE\n"
    else
        printf "Read 2 file is not AVAILABLE\n"
        exit 1
    fi

    if [[ \${fastq_list_file_status} == "AVAILABLE" ]]; then
        printf "FASTQ list file is AVAILABLE\n"
    else
        printf "FASTQ list file is not AVAILABLE\n"
        exit 1
    fi

    if [[ \${reference_file_status} == "AVAILABLE" ]]; then
        printf "Reference file is AVAILABLE\n"
    else
        printf "Reference file is not AVAILABLE\n"
        exit 1
    fi

    printf "readyForAnalysis:true\n" >> ${dataFile}
    """
}

process startAnalysis {
    debug true
    
    input:
    path(dataFile)

    output:
    path "data.txt", emit: dataFile

    script:
    """
    #!/bin/bash

    sample_id=\$(cat ${dataFile} | grep -o 'sampleId:.*' | cut -f2- -d:)
    read_1_file_id=\$(cat ${dataFile} | grep -o 'read1:.*' | cut -f2- -d:)
    read_2_file_id=\$(cat ${dataFile} | grep -o 'read2:.*' | cut -f2- -d:)
    fastq_list_file_id=\$(cat ${dataFile} | grep -o 'fastq_list:.*' | cut -f2- -d:)

    read_1_analysis_code=\$(cat ${dataFile} | grep -E "read1")
    read_2_analysis_code=\$(cat ${dataFile} | grep -E "read2")
    reference_analysis_code=\$(cat ${dataFile} | grep -E "ref_tar")

    user_reference=${userReference}-\${sample_id}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Starting Nextflow analysis...\n"

    analysis_response=\$(icav2 projectpipelines start nextflow ${pipelineId} \
        --user-reference \${user_reference} \
        --project-id ${projectId} \
        --storage-size ${storageSize} \
        --input \${reference_analysis_code} \
        --input ${fastqsAnalysisDataCode}:"\${read_1_file_id},\${read_2_file_id}" \
        --input ${fastqListDataCode}:\${fastq_list_file_id} \
        --parameters enable_map_align:true \
        --parameters enable_map_align_output:true \
        --parameters output_format:CRAM \
        --parameters enable_duplicate_marking:true \
        --parameters enable_variant_caller:true \
        --parameters vc_emit_ref_confidence:BP_RESOLUTION \
        --parameters vc_enable_vcf_output:true \
        --parameters enable_cnv:false \
        --parameters enable_sv:false \
        --parameters repeat_genotype_enable:false \
        --parameters enable_hla:false \
        --parameters enable_variant_annotation:false \
        --parameters output_file_prefix:"\${sample_id}")

    analysis_response_file="analysis_response.txt"
    touch \${analysis_response_file}
    echo "\${analysis_response}" > \${analysis_response_file}

    analysis_id=\$(cat \${analysis_response_file} | jq -r ".id")
    analysis_ref=\$(cat \${analysis_response_file} | jq -r ".reference")

    printf "[\${time_stamp}]: "
    printf "Writing id of analysis '\${analysis_ref}' to existing data file...\n"
    printf "analysisId:\${analysis_id}\n" >> ${dataFile}

    printf "Writing reference of analysis '\${analysis_ref}' to existing data file...\n"
    printf "analysisRef:\${analysis_ref}\n" >> ${dataFile}
    """
}

process checkAnalysisStatus {
    debug true
    
    input:
    path(dataFile)
    val(analysisStatusCheckInterval)

    output:
    path "data.txt", emit: dataFile

    script:
    """
    #!/bin/bash

    analysis_status_check_count=0
    analysis_status="REQUESTED"

    analysis_id=\$(cat ${dataFile} | grep -o 'analysisId:.*' | cut -f2- -d:)
    analysis_ref=\$(cat ${dataFile} | grep -o 'analysisRef:.*' | cut -f2- -d:)
    
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Checking status of analysis with id '\${analysis_id}' every ${analysisStatusCheckInterval} seconds, until status is 'SUCCEEDED'...\n"
    while true;
    do
        ((analysis_status_check_count +=1 ))
        updated_analysis_response=\$(icav2 projectanalyses get \${analysis_id})

        printf "Checking status of analysis with reference '\${analysis_ref}'...\n"
        analysis_status=\$(echo \${updated_analysis_response} | jq -r ".status")

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Current status of analysis is '\${analysis_status}'...\n"

        if [[ \${analysis_status} == "SUCCEEDED" ]]; then
            printf "Analysis SUCCEEDED\n"
            printf "analysisStatus:SUCCEEDED\n" >> ${dataFile}
            break;

        elif [[ \${analysis_status} == "FAILED" ]]; then
            printf "Analysis FAILED \n"
            printf "analysisStatus:FAILED\n" >> ${dataFile}
            exit 1

        elif [[ \${analysis_status} == "FAILED_FINAL" ]]; then
            printf "Analysis FAILED_FINAL\n"
            printf "analysisStatus:FAILED_FINAL\n" >> ${dataFile}
            exit 1

        elif [[ \${analysis_status} == "ABORTED" ]]; then
            printf "Analysis ABORTED\n"
            printf "analysisStatus:ABORTED\n" >> ${dataFile}
            exit 1

        elif [[ \${analysis_status_check_count} -gt ${analysisStatusCheckLimit} ]]; then
            printf "Analysis status has been checked more than ${analysisStatusCheckLimit} times. Stopping...\n"
            printf "analysisStatus:TIMEOUT\n" >> ${dataFile}
            exit 1

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
    path(dataFile)
    val(localDownloadPath)

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
    sample_id=\$(cat ${dataFile} | grep -o 'sampleId:.*' | cut -f2- -d:)
    read_1_file_id=\$(cat ${dataFile} | grep -o 'read1:.*' | cut -f2- -d:)
    read_2_file_id=\$(cat ${dataFile} | grep -o 'read2:.*' | cut -f2- -d:)
    analysis_output_folder_id=\$(cat ${dataFile} | grep -o 'outputFolderId:.*' | cut -f2- -d:)

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting uploaded read 1 file with ID '\${read_1_file_id}'...\n"
    icav2 projectdata delete \${read_1_file_id}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting uploaded read 2 file with ID '\${read_2_file_id}'...\n"
    icav2 projectdata delete \${read_2_file_id}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting analysis output folder with ID '\${analysis_output_folder_id}'...\n"
    icav2 projectdata delete \${analysis_output_folder_id}

    ica_upload_path="/fastq/\${sample_id}/"
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Deleting folder containing pair of FASTQ files with path '\${ica_upload_path}'...\n"
    icav2 projectdata delete \${ica_upload_path} --project-id ${projectId}

    printf "Uploaded files and analysis output folder successfully deleted.\n"
    """
}

workflow {
    fastqFilePairs = Channel.fromFilePairs(readsPairFilesUploadPath, checkIfExists:true)

    uploadFastqFilePairs(fastqFilePairs, params.projectId)
    uploadFastqFileList(uploadFastqFilePairs.out.dataFile)
    getReferenceFile(uploadFastqFileList.out.dataFile)
    checkFileStatus(getReferenceFile.out.dataFile)
    startAnalysis(checkFileStatus.out.dataFile)
    checkAnalysisStatus(startAnalysis.out.dataFile, params.analysisStatusCheckInterval)
    downloadAnalysisOutput(checkAnalysisStatus.out.dataFile, params.localDownloadPath)
    deleteData(downloadAnalysisOutput.out.dataFile)
}
