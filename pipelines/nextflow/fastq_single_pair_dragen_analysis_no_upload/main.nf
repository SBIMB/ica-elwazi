#!/usr/bin/env nextflow
nextflow.enable.dsl=2
projectId = params.projectId
read1AnalysisDataCode = params.read1AnalysisDataCode
read2AnalysisDataCode = params.read2AnalysisDataCode
fastqListDataCode = params.fastqListDataCode
referenceAnalysisDataCode = params.referenceAnalysisDataCode
pipelineId = params.pipelineId
pipelineCode = params.pipelineCode
userReference = params.userReference
storageSize = params.storageSize
referenceDirectory = params.referenceDirectory
intermediateResultsDirectory = params.intermediateResultsDirectory
fileUploadStatusCheckInterval = params.fileUploadStatusCheckInterval
analysisStatusCheckInterval = params.analysisStatusCheckInterval
analysisStatusCheckLimit = params.analysisStatusCheckLimit
sampleId = params.sampleId
read1FileId = params.read1FileId
read2FileId = params.read2FileId
fastqListFileId = params.fastqListFileId
fastqsDataCode = params.fastqsDataCode
referenceFileId = params.referenceFileId
localDownloadPath = params.localDownloadPath

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

    fastq_list_file_data_response=\$(icav2 projectdata get ${fastqListFileId})
    fastq_list_file_status=\$(echo \${fastq_list_file_data_response} | jq -r ".details.status")

    reference_file_data_response=\$(icav2 projectdata get ${referenceFileId})
    reference_file_status=\$(echo \${reference_file_data_response} | jq -r ".details.status")

    if [[ \${read_1_file_status} == "AVAILABLE" ]]; then
        printf "Read 1 file is AVAILABLE\n"
        read_1_analysis_code="${read1AnalysisDataCode}:${read1FileId}"
    else
        printf "Read 1 file is not AVAILABLE\n"
        exit 1
    fi

    if [[ \${read_2_file_status} == "AVAILABLE" ]]; then
        printf "Read 2 file is AVAILABLE\n"
        read_2_analysis_code="${read2AnalysisDataCode}:${read2FileId}"
    else
        printf "Read 2 file is not AVAILABLE\n"
        exit 1
    fi

    if [[ \${fastq_list_file_status} == "AVAILABLE" ]]; then
        printf "FASTQ list file is AVAILABLE\n"
        fastq_list_analysis_code="${fastqListDataCode}:${fastqListFileId}"
    else
        printf "FASTQ list file is not AVAILABLE\n"
        exit 1
    fi

    if [[ \${reference_file_status} == "AVAILABLE" ]]; then
        printf "Reference file is AVAILABLE\n"
        reference_file_analysis_code="${referenceAnalysisDataCode}:${referenceFileId}"
    else
        printf "Reference file is not AVAILABLE\n"
        exit 1
    fi

    data_file="data.txt"

    if ! [ -f \${data_file} ]; then
        echo "Data file does not exist. Creating one..."
        touch \${data_file}
    fi

    printf "[\${time_stamp}]: "
    printf "Writing file data to existing data file...\n"

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
    path "data.txt", emit: dataFile

    script:
    """
    #!/bin/bash

    sample_id=\$(cat ${dataFile} | grep -o 'sampleId:.*' | cut -f2- -d:)
    read_1_file_id=\$(cat ${dataFile} | grep -o 'read1:.*' | cut -f2- -d:)
    read_2_file_id=\$(cat ${dataFile} | grep -o 'read2:.*' | cut -f2- -d:)

    read_1_analysis_code=\$(cat ${dataFile} | grep -E "read1")
    read_2_analysis_code=\$(cat ${dataFile} | grep -E "read2")
    reference_analysis_code=\$(cat ${dataFile} | grep -E "ref_tar")

    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Starting Nextflow analysis...\n"
    analysis_response=\$(icav2 projectpipelines start nextflow ${pipelineId} \
        --user-reference ${userReference} \
        --project-id ${projectId} \
        --storage-size ${storageSize} \
        --input \${reference_analysis_code} \
        --input ${fastqsDataCode}:"\${read_1_file_id},\${read_2_file_id}" \
        --input ${fastqListDataCode}:${fastqListFileId} \
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
        --parameters output_file_prefix:"${sampleId}")

    analysis_response_file="analysis_response.txt"
    touch \${analysis_response_file}
    echo "\${analysis_response}" > \${analysis_response_file}

    analysis_id=\$(cat \${analysis_response_file} | jq -r ".id")
    analysis_ref=\$(cat \${analysis_response_file} | jq -r ".reference")

    printf "[\${time_stamp}]: "
    printf "Writing id of analysis '\${analysis_ref}' to existing data file...\n"
    printf "analysisId:\${analysis_id}\n" >> ${dataFile}

    printf "Writing reference of analysis '\${analysis_ref}' to existing data file...\n"
    printf "analysisRef:\${analysis_ref}\n" >> ${dataFile}    """
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
    checkFileStatus()
    startAnalysis(checkFileStatus.out.dataFile)
    checkAnalysisStatus(startAnalysis.out.dataFile, params.analysisStatusCheckInterval)
    downloadAnalysisOutput(checkAnalysisStatus.out.dataFile, params.localDownloadPath)
    deleteData(downloadAnalysisOutput.out.dataFile)
}
