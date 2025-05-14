#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process uploadCramFiles {
  debug true
  maxForks 4
  tag "$sampleId"

  input:
  tuple val(sampleId), file(cram)

  output:
  path "data.txt", emit: dataFile

  script:
  def cram_file = sampleId + ".cram"
  def crai_file = sampleId + ".cram.crai"
  def projectId = params.projectId
  """
  #!/bin/bash
    time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")
    cram_file_id=""
    crai_file_id=""

    cram_file_response="cram_file_response.txt"
    crai_file_response="crai_file_response.txt"

    touch \${cram_file_response}
    touch \${crai_file_response}

    printf "[\${time_stamp}]: "
    printf "Uploading .cram file '${cram_file}'... \n"
    cram_file_upload_response=\$(icav2 projectdata upload ${cram_file} \${ica_upload_path} --project-id ${projectId})
    echo "\${cram_file_upload_response}" > \${cram_file_response}

    printf "[\${time_stamp}]: "
    printf "Uploading .cram.crai file '${crai_file}'... \n"
    crai_file_upload_response=\$(icav2 projectdata upload ${crai_file} \${ica_upload_path} --project-id ${projectId})
    echo "\${crai_file_upload_response}" > \${crai_file_response}

    # id of file starts with 'fil.'
    cram_file_id=\$(cat \${cram_file_response} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')
    crai_file_id=\$(cat \${crai_file_response} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')

    cram_uploaded_file_data_response=\$(icav2 projectdata get \${cram_file_id})
    if [[ \$? != 0 ]]; then
        printf "Failed to fetch data about file with id '\${cram_file_id}'. \n"
        exit 1
    else
        cram_uploaded_file_name=\$(echo \${cram_uploaded_file_data_response} | jq -r ".details.name")
        printf "Name of uploaded file is '\${cram_uploaded_file_name}'. \n"
        cram_uploaded_file_path=\$(echo \${cram_uploaded_file_data_response} | jq -r ".details.path")
        printf "Path of uploaded .cram file is '\${cram_uploaded_file_path}'. \n"
    fi

    crai_uploaded_file_data_response=\$(icav2 projectdata get \${crai_file_id})
    if [[ \$? != 0 ]]; then
        printf "Failed to fetch data about file with id '\${crai_file_id}'. \n"
        exit 1
    else
        crai_uploaded_file_name=\$(echo \${crai_uploaded_file_data_response} | jq -r ".details.name")
        printf "Name of uploaded file is '\${crai_uploaded_file_name}'. \n"
        crai_uploaded_file_path=\$(echo \${crai_uploaded_file_data_response} | jq -r ".details.path")
        printf "Path of uploaded .crai file is '\${crai_uploaded_file_path}'. \n"
    fi

    data_file="data.txt"

    if ! [ -f \${data_file} ]; then
        echo "Data file does not exist. Creating one..."
        touch \${data_file}
    fi

    printf "[\${time_stamp}]: "
    printf "Writing file data to existing data file...\n"

    printf "sampleId:${sampleId}\n" >> \${data_file}
    printf "cram:\${cram_file_id}\n" >> \${data_file}
    printf "cramIndex:\${crai_file_id}\n" >> \${data_file}
  """
}

process getReferenceFile {
  debug true

  input:
  path(dataFile)

  output:
  path "data.txt", emit: dataFile

  script:
  def projectId = params.projectId
  def referenceAnalysisDataCode = params.referenceAnalysisDataCode
  def referenceFileId = params.referenceFileId
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
    def fileStatusCheckInterval = params.fileStatusCheckInterval
    def fileStatusCheckLimit = params.fileStatusCheckLimit
    """
    #!/bin/bash
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Getting status of uploaded files...\n"

    cram_file_id=\$(cat ${dataFile} | grep -o 'cram:.*' | cut -f2- -d:)
    crai_file_id=\$(cat ${dataFile} | grep -o 'cramIndex:.*' | cut -f2- -d:)
    fastq_list_file_id=\$(cat ${dataFile} | grep -o 'fastq_list:.*' | cut -f2- -d:)
    reference_file_id=\$(cat ${dataFile} | grep -o 'ref_tar:.*' | cut -f2- -d:)

    file_status_check_count=0

    while true;
    do
        ((file_status_check_count +=1 ))

        cram_file_data_response=\$(icav2 projectdata get \${cram_file_id})

        printf "Checking status of cram file with id '\${cram_file_id}'...\n"
        cram_file_status=\$(echo \${cram_file_data_response} | jq -r ".details.status")
        printf "[\${timeStamp}]: Current status of cram file is '\${cram_file_status}'...\n"

        crai_file_data_response=\$(icav2 projectdata get \${crai_file_id})

        printf "Checking status of crai file with id '\${crai_file_id}'...\n"
        crai_file_status=\$(echo \${crai_file_data_response} | jq -r ".details.status")

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Current status of crai file is '\${crai_file_status}'...\n"

        fastq_list_file_data_response=\$(icav2 projectdata get \${fastq_list_file_id})

        printf "Checking status of fastq list file with id '\${fastq_list_file_id}'...\n"
        fastq_list_file_status=\$(echo \${fastq_list_file_data_response} | jq -r ".details.status")

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Current status of fastq list file is '\${fastq_list_file_status}'...\n"

        if [ \${cram_file_status} == "AVAILABLE" ] && [ \${crai_file_status} == "AVAILABLE" ] && [ \${fastq_list_file_status} == "AVAILABLE" ]; then
            printf "CRAM file is AVAILABLE\n"
            printf "CRAM index file is AVAILABLE\n"
            printf "FASTQ list file is AVAILABLE\n"
            break;

        elif [ \${file_status_check_count} -gt ${fileStatusCheckLimit} ]; then
            printf "File status has been checked more than ${fileStatusCheckLimit} times. Stopping...\n"
            exit 1

        else
            printf "File availability still in progress...\n"
        fi

        sleep ${fileStatusCheckInterval};
    done

    reference_file_data_response=\$(icav2 projectdata get \${reference_file_id})
    reference_file_status=\$(echo \${reference_file_data_response} | jq -r ".details.status")

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
    def projectId = params.projectId
    def cramAnalysisDataCode = params.cramAnalysisDataCode
    def cramIndexAnalysisDataCode = params.cramIndexAnalysisDataCode
    def pipelineId = params.pipelineId
    def userReference = params.userReference
    def storageSize = params.storageSize
    """
    #!/bin/bash

    sample_id=\$(cat ${dataFile} | grep -o 'sampleId:.*' | cut -f2- -d:)
    cram_file_id=\$(cat ${dataFile} | grep -o 'cram:.*' | cut -f2- -d:)
    crai_file_id=\$(cat ${dataFile} | grep -o 'cramIndex:.*' | cut -f2- -d:)
    fastq_list_file_id=\$(cat ${dataFile} | grep -o 'fastq_list:.*' | cut -f2- -d:)

    cram_analysis_code=\$(cat ${dataFile} | grep -E "cram")
    crai_analysis_code=\$(cat ${dataFile} | grep -E "cramIndex")
    reference_analysis_code=\$(cat ${dataFile} | grep -E "ref_tar")

    user_reference=${userReference}-\${sample_id}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Starting Nextflow analysis...\n"

    analysis_response=\$(icav2 projectpipelines start nextflow ${pipelineId} \
        --user-reference \${user_reference} \
        --project-id ${projectId} \
        --storage-size ${storageSize} \
        --input \${reference_analysis_code} \
        --input ${cramAnalysisDataCode}:"\${cram_file_id}\
        --input ${cramIndexAnalysisDataCode}:\${fastq_list_file_id} \
        --parameters enable_map_align:true \
        --parameters enable_map_align_output:true \
        --parameters output_format:CRAM \
        --parameters enable_duplicate_marking:true \
        --parameters enable_variant_caller:true \
        --parameters vc_emit_ref_confidence:GVCF \
        --parameters vc_enable_vcf_output:true \
        --parameters enable_cnv:true \
        --parameters enable_sv:true \
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

    output:
    path "data.txt", emit: dataFile

    script:
    def analysisStatusCheckInterval = params.analysisStatusCheckInterval
    def analysisStatusCheckLimit = params.analysisStatusCheckLimit
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

    output:
    path "data.txt", emit: dataFile

    script:
    def localDownloadPath = params.localDownloadPath
    """
    #!/bin/bash
    analysis_status=\$(cat ${dataFile} | grep -o 'analysisStatus:.*' | cut -f2- -d:)
    analysis_id=\$(cat ${dataFile} | grep -o 'analysisId:.*' | cut -f2- -d:)

    if [ "\$analysis_status" != "SUCCEEDED" ]; then
        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Analysis did not succeed in previous process. Data will not be downloaded nor deleted...\n"
        printf "deleteData:false\n" >> ${dataFile}
    else
        printf "[\${time_stamp}]: "
        printf "Fetching analysis output response...\n"
        analysis_output_response=\$(icav2 projectanalyses output \${analysis_id})
        analysis_output_folder_id=\$(echo \${analysis_output_response} | jq -r ".items[].data[].dataId")
        printf "Analysis output folder ID is '\${analysis_output_folder_id}'\n"
        printf "Writing id of analysis output folder to existing data file...\n"
        printf "outputFolderId:\${analysis_output_folder_id}\n" >> ${dataFile}
        printf "deleteData:true\n" >> ${dataFile}

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Downloading analysis output folder with ID '\${analysis_output_folder_id}' to '${localDownloadPath}'...\n"

        icav2 projectdata download \${analysis_output_folder_id} ${localDownloadPath}
    fi
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
    delete_data=\$(cat ${dataFile} | grep -o 'deleteData:.*' | cut -f2- -d:)

    if [ "\$delete_data" = "false" ]; then
        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Data will NOT be deleted due to failed analysis...\n"
    else
        sample_id=\$(cat ${dataFile} | grep -o 'sampleId:.*' | cut -f2- -d:)
        cram_file_id=\$(cat ${dataFile} | grep -o 'cram:.*' | cut -f2- -d:)
        crai_file_id=\$(cat ${dataFile} | grep -o 'cramIndex:.*' | cut -f2- -d:)
        fastq_list_file_id=\$(cat ${dataFile} | grep -o 'fastq_list:.*' | cut -f2- -d:)
        analysis_output_folder_id=\$(cat ${dataFile} | grep -o 'outputFolderId:.*' | cut -f2- -d:)

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Deleting uploaded cram file with ID '\${cram_file_id}'...\n"
        icav2 projectdata delete \${cram_file_id}

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Deleting uploaded crai file with ID '\${crai_file_id}'...\n"
        icav2 projectdata delete \${crai_file_id}

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Deleting uploaded CSV file with ID '\${fastq_list_file_id}'...\n"
        icav2 projectdata delete \${fastq_list_file_id}

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Deleting analysis output folder with ID '\${analysis_output_folder_id}'...\n"
        icav2 projectdata delete \${analysis_output_folder_id}

        printf "Uploaded files and analysis output folder successfully deleted.\n"
    fi
    """
}

workflow {
    cramFilePairsChannel = Channel.fromFilePairs(params.cramFilePairsUploadPath, checkIfExists:true) { 
      file -> file.name.replaceAll(/.cram|.crai$/,'') 
    }
    uploadCramFiles(cramFilePairsChannel)
    getReferenceFile(uploadCramFiles.out.dataFile)
    checkFileStatus(getReferenceFile.out.dataFile)
    // startAnalysis(checkFileStatus.out.dataFile)
    // checkAnalysisStatus(startAnalysis.out.dataFile)
    // downloadAnalysisOutput(checkAnalysisStatus.out.dataFile)
    // deleteData(downloadAnalysisOutput.out.dataFile)
}
