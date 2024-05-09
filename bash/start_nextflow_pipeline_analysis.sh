#!/bin/bash

project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
pipeline_code="basic_pipeline"
pipeline_id="bfecca03-6443-45bd-b313-e4f555cd0748"
user_reference="regan_test_analysis_03"
storage_size="Small"
file_name="NZ_GG704940.fa"
file_id="fil.09b47f75e6014039ac4908dc6dfe2358"

printf "Fetching list of project analyses... \n"
projectanalyses_list_response=$(icav2 projectanalyses list --project-id $project_id)

if [[ $? != 0 ]]; then
    printf "Failed to fetch projectanalyses list. \n"
    exit 1
else
    analysis_id=$(echo $projectanalyses_response | jq -r ".items[] | select(.userReference == \"$user_reference\").id")
    printf "Fetching project analysis with id '$analysis_id'... \n"
    analysis_response=$(icav2 projectanalyses input $analysis_id)
    if [[ $? != 0 ]]; then
        printf "Failed to projectanalysis with id '$analysis_id'. \n"
        exit 1
    else
        analysis_data_id=$(echo $analysis_response | jq -r ".items[].analysisData[] | select(.name == \"$file_name\").dataId")
        analysis_data_code=$(echo $analysis_response | jq -r ".items[].code")

        printf "Constructing 'analysisCode:dataId' from analysis response... \n"
        analysis_ref="$analysis_data_code:$analysis_data_id"

        printf "Constructing 'analysisCode:fileId' from analysis response and file details... \n"
        file_ref="$analysis_data_code:$file_id"

        printf "Starting Nextflow analysis... \n"
        icav2 projectpipelines start nextflow $pipeline_id \
            --user-reference $user_reference \
            --project-id $project_id \
            --storage-size $storage_size \
            --input $file_ref
    fi
fi






