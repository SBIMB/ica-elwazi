#! /bin/bash

project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
analysis_user_reference="regan_test_analysis_03"

analysis_status="REQUESTED"
interval_in_seconds=240

printf "Polling ICA V2 every $interval_in_seconds seconds, until status of analysis is 'SUCCEEDED'\n"

while true;
do
    printf "Fetching list of project analyses... \n"
    projectanalyses_response=$(icav2 projectanalyses list --project-id $project_id)
    if [[ $? == 0 ]]; then
        printf "Fetching status of project analysis with user reference '$analysis_user_reference'... \n"
        analysis_status=$(echo $projectanalyses_response | jq -r ".items[] | select(.userReference == \"$analysis_user_reference\").status")
        time_stamp=$(date +"%Y-%m-%d_%H-%M-%S")
        printf "$time_stamp: "
        if [ $analysis_status == "SUCCEEDED" ]; then
            printf "Analysis SUCCEEDED \n"
            break;

        elif [ $analysis_status == "FAILED" ]; then
            printf "Analysis FAILED \n"
            break;

        elif [ $analysis_status == "FAILED_FINAL" ]; then
            printf "Analysis FAILED_FINAL \n"
            break;

        elif [ $analysis_status == "ABORTED" ]; then
            printf "Analysis ABORTED \n"
            break;
        else
            printf "Analysis still in progress... \n"
        fi
    else
        printf "Failed to fetch status of project analysis with user reference '$analysis_user_reference'... \n"
        break;
    fi
    sleep $interval_in_seconds;
done