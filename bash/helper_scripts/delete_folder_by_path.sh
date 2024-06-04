#!/bin/bash

folder_id="fol.eaea710fdc0b448d2f8608dc6c0c20cf"
folder_path="/Citrobacter_freundii_4_7_47CFAA_uid46379/"
folder_name="Citrobacter_freundii_4_7_47CFAA_uid46379"
local_storage_path="$HOME/Documents/ica_data_downloads/"
project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
time_stamp=$(date +"%Y-%m-%d_%H-%M-%S")

printf "Checking if folder exists in ICA storage... \n"
projectdata_list_response=$(icav2 projectdata list --project-id $project_id)

printf "$time_stamp: "
if [[ $? != 0 ]]; then
    printf "Failed to fetch projectdata list. \n"
else
    folder_ica_storage_path=$(echo $projectdata_list_response | jq -r ".items[].details | select(.name == \"$folder_name\").path")
    if [[ $folder_ica_storage_path == "" ]]; then
        printf "Folder with name '$folder_name' does NOT exist in ICA storage. \n"
        exit 1
    else
        printf "Folder with name '$folder_name' exists in ICA storage at '$folder_ica_storage_path'. \n"
        printf "Checking if folder exists in local storage... \n"
        folder_local_storage_path="${local_storage_path}/${folder_name}"
        if [[ -d "$folder_local_storage_path" ]]; then
            printf "Folder '$folder_name' exists in local path '$local_storage_path'. Deleting from ICA storage... \n"
            icav2 projectdata delete $folder_ica_storage_path --project-id $project_id
            if [[ $? != 0 ]]; then
                printf "Failed to delete folder. \n"
            else
                printf "Folder '$folder_name' successfully deleted from ICA storage path '$folder_ica_storage_path'. \n"
            fi
        else
            printf "Folder '$folder_name' does NOT exist in local path '$local_storage_path'. Delete process stopped... \n"
            exit 1
        fi
    fi
fi

