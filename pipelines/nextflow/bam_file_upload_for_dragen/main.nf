#!/usr/bin/env nextflow
nextflow.enable.dsl=2
projectId = params.projectId
bamAnalysisDataCode = params.bamAnalysisDataCode
referenceAnalysisDataCode = params.referenceAnalysisDataCode
pipelineId = params.pipelineId
pipelineCode = params.pipelineCode
userReference = params.userReference
storageSize = params.storageSize
fileUploadStatusCheckInterval = params.fileUploadStatusCheckInterval
analysisStatusCheckInterval = params.analysisStatusCheckInterval
bamFilesUploadPath = params.bamFilesUploadPath
bamFilePairsUploadPath = params.bamFilePairsUploadPath
referenceFileUploadPath = params.referenceFileUploadPath
localDownloadPath = params.localDownloadPath
icaUploadPath = params.icaUploadPath

process uploadBamFiles {
  debug true
  maxForks 3
  tag "$sampleId"

  input:
  tuple val(sampleId), file(bam)

  output:
  path "manifest.tsv", emit: manifestFile

  script:
  def bam_file = sampleId + ".bam"
  def bai_file = sampleId + ".bam.bai"
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

  manifest_file="manifest.tsv"
  manifest_tab=\$(printf "\t")

  if ! [ -f \${manifest_file} ]; then
      echo "Manifest file does not exist. Creating one..."
      touch \${manifest_file}

      manifest_header_1="file_name"
      manifest_header_2="file_analysis_code"

      printf "[\${time_stamp}]: "
      printf "Writing file data to manifest...\n"

      printf "\${manifest_header_1} \${manifest_tab} \${manifest_header_2}\n" >> \${manifest_file}
      printf "\${bam_file_id} \${manifest_tab} ${bamAnalysisDataCode}:\${bam_file_id}\n" >> \${manifest_file}
  else
      printf "[\${time_stamp}]: "
      printf "Writing file data to existing manifest...\n"

      printf "\${bam_file_id} \${manifest_tab} ${bamAnalysisDataCode}:\${bam_file_id}\n" >> \${manifest_file}
  fi
  """
}

workflow {
  bamFilePairsChannel = Channel.fromFilePairs(params.bamFilePairsUploadPath, checkIfExists:true) { 
    file -> file.name.replaceAll(/.bam|.bai$/,'') 
  }
  uploadBamFiles(bamFilePairsChannel)
}
