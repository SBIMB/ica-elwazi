#!/bin/bash

project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
pipeline_id="858708e9-6cdd-42ae-ad97-8b86ee07cdc5"
analysis_id="f046c8e2-4add-4716-acd3-c5674b1fed74"
user_reference="regan_dragen_germline_whole_genome_test"
analysis_status="REQUESTED"
analysis_status_check_count=0
analysis_status_check_limit=1000
interval_in_seconds=300

printf "Polling ICA V2 every $interval_in_seconds seconds, until status of analysis is 'SUCCEEDED'\n"

timeStamp=$(date +"%Y-%m-%d %H:%M:%S")
printf "[$timeStamp]: Checking status of analysis with id '$analysis_id' every $interval_in_seconds seconds, until status is 'SUCCEEDED'...\n"
while true;
do
    printf "Checking the analysis status: '$analysis_status_check_count' \n"
    ((analysis_status_check_count += 1))
    updated_analysis_response=$(icav2 projectanalyses get $analysis_id --project-id $project_id)

    printf "Checking status of analysis with user reference '$user_reference'...\n"
    analysis_status=$(echo $updated_analysis_response | jq -r ".status")

    timeStamp=$(date +"%Y-%m-%d %H:%M:%S")
    printf "[$timeStamp]: Current status of analysis is '$analysis_status'...\n"

    if [[ $analysis_status == "SUCCEEDED" ]]; then
        printf "Analysis SUCCEEDED\n"
        printf "Fetching analysis output response...\n"
        analysis_output_response=$(icav2 projectanalyses output $analysis_id)
        analysis_output_folder_id=$(echo $analysis_output_response | jq -r ".items[].data[].dataId")
        printf "Analysis output folder ID is '$analysis_output_folder_id'\n"
        break;

    elif [[ $analysis_status == "FAILED" ]]; then
        printf "Analysis FAILED \n"
        break;

    elif [[ $analysis_status == "FAILED_FINAL" ]]; then
        printf "Analysis FAILED_FINAL\n"
        break;

    elif [[ $analysis_status == "ABORTED" ]]; then
        printf "Analysis ABORTED\n"
        break;

    elif [[ $analysis_status_check_count -gt $analysis_status_check_limit ]]; then
        printf "Analysis status has been checked more than $analysis_status_check_limit times. Stopping...\n"
        break;

    else
        printf "Analysis still in progress...\n"
    fi

    sleep $interval_in_seconds;
done