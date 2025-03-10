#!/bin/bash

analysis_status_check_count=0
analysis_status="REQUESTED"

analysis_id=""
analysis_ref=""

stop_polling=false

list_of_failed_analysis_status_values=("FAILED", "FAILED_FINAL", "ABORTED")

time_stamp=$(date +"%Y-%m-%d %H:%M:%S")
printf "[$time_stamp]: Checking status of analysis with id '$analysis_id' every $analysis_status_check_interval seconds, until status is 'SUCCEEDED'...\n"
while true;
do
    ((analysis_status_check_count +=1 ))
    updated_analysis_response=$(icav2 projectanalyses get $analysis_id)

    printf "Checking status of analysis with reference '$analysis_ref'...\n"
    analysis_status=$(echo $updated_analysis_response | jq -r ".status")

    time_stamp=$(date +"%Y-%m-%d %H:%M:%S")
    printf "[$time_stamp]: Current status of analysis is '$analysis_status'...\n"

    if [[ $analysis_status == "SUCCEEDED" ]]; then
        printf "Analysis SUCCEEDED\n"
        stop_polling=true

    elif [[ $analysis_status_check_count -gt $analysis_status_check_limit ]]; then
        printf "Analysis status has been checked more than $analysis_status_check_limit times. Stopping...\n"
        stop_polling=true

    else
        printf "Analysis still in progress...\n"
        stop_polling=false
    fi

    for status in $list_of_failed_analysis_status_values
    do
        if [[ $analysis_status == "$status" ]]; then 
            printf "Analysis $status\n"
            stop_polling=true
        fi
    done;

    if $stop_polling; then break;

    sleep $analysis_status_check_interval;
done

printf "[$time_stamp]: Status of analysis with reference '$analysis_ref' is $analysis_status.\n"
