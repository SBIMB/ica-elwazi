#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process uploadBamFiles {
  debug true
  maxForks 4
  tag "$sampleId"

  input:
  tuple val(sampleId), file(bamPair)

  output:
  path "data.txt", emit: dataFile

  script:
  def projectId = params.projectId
  def bamAnalysisDataCode = params.bamAnalysisDataCode
  def bamIndexAnalysisDataCode = params.bamIndexAnalysisDataCode
  def (bam_file, bai_file) = bamPair
  """
  #!/bin/bash
    time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")
    bam_file_id=""
    bai_file_id=""

    bam_file_response="bam_file_response.txt"
    bai_file_response="bai_file_response.txt"

    touch \${bam_file_response}
    touch \${bai_file_response}

    printf "[\${time_stamp}]: "
    printf "Uploading .bam file '${bam_file}'... \n"
    bam_file_upload_response=\$(icav2 projectdata upload ${bam_file} \${ica_upload_path} --project-id ${projectId})
    echo "\${bam_file_upload_response}" > \${bam_file_response}

    printf "[\${time_stamp}]: "
    printf "Uploading .bam.bai file '${bai_file}'... \n"
    bai_file_upload_response=\$(icav2 projectdata upload ${bai_file} \${ica_upload_path} --project-id ${projectId})
    echo "\${bai_file_upload_response}" > \${bai_file_response}

    # id of file starts with 'fil.'
    bam_file_id=\$(cat \${bam_file_response} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')
    bai_file_id=\$(cat \${bai_file_response} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')

    bam_uploaded_file_data_response=\$(icav2 projectdata get \${bam_file_id} --project-id ${projectId})
    bam_uploaded_file_name=\$(echo \${bam_uploaded_file_data_response} | jq -r ".details.name")
    printf "Name of uploaded file is '\${bam_uploaded_file_name}'. \n"

    bai_uploaded_file_data_response=\$(icav2 projectdata get \${bai_file_id} --project-id ${projectId})
    bai_uploaded_file_name=\$(echo \${bai_uploaded_file_data_response} | jq -r ".details.name")
    printf "Name of uploaded file is '\${bai_uploaded_file_name}'. \n"

    data_file="data.txt"

    if ! [ -f \${data_file} ]; then
        echo "Data file does not exist. Creating one..."
        touch \${data_file}
    fi

    printf "[\${time_stamp}]: "
    printf "Writing file data to existing data file...\n"

    printf "sampleId:${sampleId}\n" >> \${data_file}
    printf "${bamAnalysisDataCode}:\${bam_file_id}\n" >> \${data_file}
    printf "${bamIndexAnalysisDataCode}:\${bai_file_id}\n" >> \${data_file}
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

    bam_file_id=\$(cat ${dataFile} | grep -o 'bam:.*' | cut -f2- -d:)
    bai_file_id=\$(cat ${dataFile} | grep -o 'bamIndex:.*' | cut -f2- -d:)
    reference_file_id=\$(cat ${dataFile} | grep -o 'ref_tar:.*' | cut -f2- -d:)

    file_status_check_count=0

    while true;
    do
        ((file_status_check_count +=1 ))

        bam_file_data_response=\$(icav2 projectdata get \${bam_file_id})

        printf "Checking status of bam file with id '\${bam_file_id}'...\n"
        bam_file_status=\$(echo \${bam_file_data_response} | jq -r ".details.status")
        printf "[\${timeStamp}]: Current status of bam file is '\${bam_file_status}'...\n"

        bai_file_data_response=\$(icav2 projectdata get \${bai_file_id})

        printf "Checking status of bai file with id '\${bai_file_id}'...\n"
        bai_file_status=\$(echo \${bai_file_data_response} | jq -r ".details.status")

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Current status of bai file is '\${bai_file_status}'...\n"

        if [[ (\${bam_file_status} == "AVAILABLE") && (\${bai_file_status} == "AVAILABLE") ]]; then
            printf "bam file is AVAILABLE\n"
            printf "bam index file is AVAILABLE\n"
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
    def bamAnalysisDataCode = params.bamAnalysisDataCode
    def bamIndexAnalysisDataCode = params.bamIndexAnalysisDataCode
    def pipelineId = params.pipelineId
    def userReference = params.userReference
    def storageSize = params.storageSize
    """
    #!/bin/bash

    sample_id=\$(cat ${dataFile} | grep -o 'sampleId:.*' | cut -f2- -d:)
    bam_file_id=\$(cat ${dataFile} | grep -o 'bam:.*' | cut -f2- -d:)
    bai_file_id=\$(cat ${dataFile} | grep -o 'bamIndex:.*' | cut -f2- -d:)

    bam_analysis_code=\$(cat ${dataFile} | grep -E "bam")
    bai_analysis_code=\$(cat ${dataFile} | grep -E "bamIndex")
    reference_analysis_code=\$(cat ${dataFile} | grep -E "ref_tar")

    user_reference=${userReference}-\${sample_id}

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Starting Nextflow analysis...\n"

    analysis_response=\$(icav2 projectpipelines start nextflow ${pipelineId} \
        --user-reference \${user_reference} \
        --project-id ${projectId} \
        --storage-size ${storageSize} \
        --input \${reference_analysis_code} \
        --input ${bamAnalysisDataCode}:\${bam_file_id} \
        --input ${bamIndexAnalysisDataCode}:\${bai_file_id} \
        --parameters enable_map_align:false \
        --parameters enable_map_align_output:false \
        --parameters enable_duplicate_marking:false \
        --parameters enable_variant_caller:true \
        --parameters vc_emit_ref_confidence:GVCF \
        --parameters vc_enable_vcf_output:true \
        --parameters enable_cnv:true \
        --parameters enable_sv:true \
        --parameters repeat_genotype_enable:true \
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
        bam_file_id=\$(cat ${dataFile} | grep -o 'bam:.*' | cut -f2- -d:)
        bai_file_id=\$(cat ${dataFile} | grep -o 'bamIndex:.*' | cut -f2- -d:)
        analysis_output_folder_id=\$(cat ${dataFile} | grep -o 'outputFolderId:.*' | cut -f2- -d:)

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Deleting uploaded bam file with ID '\${bam_file_id}'...\n"
        icav2 projectdata delete \${bam_file_id}

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Deleting uploaded bai file with ID '\${bai_file_id}'...\n"
        icav2 projectdata delete \${bai_file_id}

        timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
        printf "[\${timeStamp}]: Deleting analysis output folder with ID '\${analysis_output_folder_id}'...\n"
        icav2 projectdata delete \${analysis_output_folder_id}

        printf "Uploaded files and analysis output folder successfully deleted.\n"
    fi
    """
}

workflow {
    def bamFilePairsUploadPath = params.bamFilePairsUploadPath
    bamFilePairsChannel = Channel.fromFilePairs(bamFilePairsUploadPath, checkIfExists:true) { 
        file -> file.name.replaceAll(/.bam|.bai$/,'')
    }
    uploadBamFiles(bamFilePairsChannel)
    getReferenceFile(uploadBamFiles.out.dataFile)
    checkFileStatus(getReferenceFile.out.dataFile)
    startAnalysis(checkFileStatus.out.dataFile)
    checkAnalysisStatus(startAnalysis.out.dataFile)
    downloadAnalysisOutput(checkAnalysisStatus.out.dataFile)
    deleteData(downloadAnalysisOutput.out.dataFile)
}
