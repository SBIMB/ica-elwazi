#!/usr/bin/env nextflow
nextflow.enable.dsl=2
projectId = params.projectId
bamAnalysisDataCode = params.bamAnalysisDataCode
bamIndexAnalysisDataCode = params.bamIndexAnalysisDataCode
referenceAnalysisDataCode = params.referenceAnalysisDataCode
pipelineId = params.pipelineId
pipelineCode = params.pipelineCode
userReference = params.userReference
storageSize = params.storageSize
fileUploadStatusCheckInterval = params.fileUploadStatusCheckInterval
analysisStatusCheckInterval = params.analysisStatusCheckInterval
sampleId = params.sampleId
bamFilesUploadPath = params.bamFilesUploadPath
bamFilePairsUploadPath = params.bamFilePairsUploadPath
referenceFileUploadPath = params.referenceFileUploadPath
referenceFileIcaPath = params.referenceFileIcaPath
referenceFileId = params.referenceFileId
localDownloadPath = params.localDownloadPath

process uploadBamFileIndexPair {
  debug true

  input:
  val(bamFiles)

  output:
  path "data.txt", emit: dataFile

  script:
  def bam_file = bamFiles.find { file -> file.name.endsWith(".bam") }
  def bai_file = bamFiles.find { file -> file.name.endsWith(".bam.bai") }
  def bam_file_name = bam_file.baseName + ".bam"
  def bai_file_name = bai_file.baseName + ".bai"

  """
  #!/bin/bash
  time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")
  bam_file_id=""
  bai_file_id=""
  ica_upload_path="/bam/$sampleId/"

  printf "sampleId: $sampleId \n"
  printf "bam_file: $bam_file \n"
  printf "bai_file: $bai_file \n"

  bam_file_response="bam_file_response.txt"
  bai_file_response="bai_file_response.txt"

  touch \${bam_file_response}
  touch \${bai_file_response}

  printf "[\${time_stamp}]: "
  printf "Uploading .bam file '${bam_file_name}'... \n"
  bam_file_upload_response=\$(icav2 projectdata upload ${bam_file} \${ica_upload_path} --project-id ${projectId})
  echo "\${bam_file_upload_response}" > \${bam_file_response}

  printf "[\${time_stamp}]: "
  printf "Uploading .bam.bai file '${bai_file_name}'... \n"
  bai_file_upload_response=\$(icav2 projectdata upload ${bai_file} \${ica_upload_path} --project-id ${projectId})
  echo "\${bai_file_upload_response}" > \${bai_file_response}

  # id of file starts with 'fil.'
  bam_file_id=\$(cat \${bam_file_response} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')
  bai_file_id=\$(cat \${bai_file_response} | grep -i '\"id\": \"fil' | grep -o 'fil.[^\"]*')

  bam_file_uploaded_file_data_response=\$(icav2 projectdata get \${bam_file_id})
  if [[ \$? != 0 ]]; then
      printf "Failed to fetch data about file with id '\${bam_file_id}'. \n"
      exit 1
  else
      bam_file_uploaded_file_path=\$(echo \${bam_file_uploaded_file_data_response} | jq -r ".details.path")
      printf "Path of uploaded .bam file is '\${bam_file_uploaded_file_path}'. \n"
  fi

  bai_file_uploaded_file_data_response=\$(icav2 projectdata get \${bai_file_id})
  if [[ \$? != 0 ]]; then
      printf "Failed to fetch data about file with id '\${bai_file_id}'. \n"
      exit 1
  else
      bai_file_uploaded_file_path=\$(echo \${bai_file_uploaded_file_data_response} | jq -r ".details.path")
      printf "Path of uploaded .bai file is '\${bai_file_uploaded_file_path}'. \n"
  fi

  data_file="data.txt"

  if ! [ -f \${data_file} ]; then
      echo "Data file does not exist. Creating one..."
      touch \${data_file}
  fi

  printf "[\${time_stamp}]: "
  printf "Writing file data to existing data file...\n"

  printf "sampleId:${sampleId}\n" >> \${data_file}
  printf "${bamAnalysisDataCode}:\${bam_file_id}\n" >> \${data_file}

  echo "${bamFiles}"
  echo "${bam_file}"
  echo "${bai_file}"

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
    bam_file_id=\$(cat ${dataFile} | grep -o 'bams:.*' | cut -f2- -d:)
    bai_file_id=\$(cat ${dataFile} | grep -o 'bais:.*' | cut -f2- -d:)

    reference_analysis_code=\$(cat ${dataFile} | grep -E "ref_tar")

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Starting Nextflow analysis...\n"

    analysis_response=\$(icav2 projectpipelines start nextflow ${pipelineId} \
        --user-reference ${userReference} \
        --project-id ${projectId} \
        --storage-size ${storageSize} \
        --input \${reference_analysis_code} \
        --input ${bamAnalysisDataCode}:"\${bam_file_id} \
        --input ${bamIndexAnalysisDataCode}:"\${bai_file_id} \
        --parameters enable_map_align:true \
        --parameters enable_variant_caller:true \
        --parameters vc_emit_ref_confidence:BP_RESOLUTION \
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
            break;

        elif [[ \${analysis_status} == "FAILED_FINAL" ]]; then
            printf "Analysis FAILED_FINAL\n"
            printf "analysisStatus:FAILED_FINAL\n" >> ${dataFile}
            break;

        elif [[ \${analysis_status} == "ABORTED" ]]; then
            printf "Analysis ABORTED\n"
            printf "analysisStatus:ABORTED\n" >> ${dataFile}
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

    printf "Uploaded file and analysis output folder successfully deleted.\n"
    """
}

workflow {
    bamFiles = Channel.fromPath(bamFilesUploadPath, checkIfExists:true)
    .filter( { file -> file.name.startsWith("${sampleId}") } )
    .toList()
    .view()
    uploadBamFileIndexPair(bamFiles)
    getReferenceFile(uploadBamFileIndexPair.out.dataFile)
    startAnalysis(getReferenceFile.out.dataFile)
    checkAnalysisStatus(startAnalysis.out.dataFile, params.analysisStatusCheckInterval)
    downloadAnalysisOutput(checkAnalysisStatus.out.dataFile, params.localDownloadPath)
    deleteData(downloadAnalysisOutput.out.dataFile)
}