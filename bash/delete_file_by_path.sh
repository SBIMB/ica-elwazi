#!/bin/bash

file_name="NZ_GG704942.fa"
local_storage_path="Documents/ica_data_downloads/"
project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
time_stamp=$(date +"%Y-%m-%d_%H-%M-%S")

printf "Checking if file exists in ICA storage... \n"
projectdata_list_response=$(icav2 projectdata list --project-id $project_id)

printf "$time_stamp: "
if [[ $? != 0 ]]; then
    printf "Failed to fetch projectdata list. \n"
else
    file_ica_storage_path=$(echo $projectdata_list_response | jq -r ".items[].details | select(.name == \"$file_name\").path")
    if [[ $file_ica_storage_path == "" ]]; then
        printf "File with name '$file_name' does NOT exist in ICA storage. \n"
        exit 1
    else
        printf "File with name '$file_name' exists in ICA storage at '$file_ica_storage_path'. \n"
        printf "Checking if file exists in local storage... \n"
        file_local_storage_path="${local_storage_path}/${file_name}"
        if [[ -f "$file_local_storage_path" ]]; then
            printf "File '$file_name' exists in local path '$local_storage_path'. Deleting from ICA storage... \n"
            icav2 projectdata delete $file_ica_storage_path --project-id $project_id
            if [[ $? != 0 ]]; then
                printf "Failed to delete file. \n"
            else
                printf "File '$file_name' successfully deleted from ICA storage path '$file_ica_storage_path'. \n"
            fi
        else
            printf "File '$file_name' does NOT exist in local path '$local_storage_path'. Delete process stopped... \n"
            exit 1
        fi
    fi
fi

