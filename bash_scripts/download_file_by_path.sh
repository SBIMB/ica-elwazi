#!/bin/bash

file_name="NZ_GG704940.fa"
local_download_path="Documents/ica_data_downloads/"
project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
interval_in_seconds=30
download_check_count=0
download_check_limit=10
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
    else
        printf "File with name '$file_name' exists in ICA storage at '$file_ica_storage_path'. \n"
        printf "Downloading file '$file_name' from ICA storage at '$file_ica_storage_path'... \n"
        file_download_response=$(icav2 projectdata download $file_ica_storage_path $local_download_path)
        if [[ $? != 0 ]]; then
            printf "Failed to send download command. \n"
        else
            file_download_path="${local_download_path}/${file_name}"
            while true;
            do
                printf "$time_stamp: "
                printf "Checking if download is complete... \n"
                ((download_check_count+=1))
                if [[ -f "$file_download_path" ]]; then
                    printf "File '$file_name' successfully saved in local path '$local_download_path'. \n"
                    break;
                elif [ $download_check_count -gt $download_check_limit ]; then
                    printf "Download check: '$download_check_count'"
                    break;
                else
                    printf "File '$file_name' is still downloading... \n"
                fi
                sleep $interval_in_seconds;
            done
        fi
    fi
fi
