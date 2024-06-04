#! /bin/bash

project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
pipeline_id="bfecca03-6443-45bd-b313-e4f555cd0748"
analysis_id="5dfa985b-3ffb-4aab-ab06-9b62514a8eb4"
local_download_path="${HOME}/Documents/ica_data_downloads/"
time_stamp=$(date +"%Y-%m-%d_%H:%M:%S")

# analysis_response=$(icav2 projectanalyses get $analysis_id)

printf "$time_stamp: "
printf "Fetching analysis output of analysis with id '$analysis_id'... \n"
analysis_output=$(icav2 projectanalyses output $analysis_id --project-id $project_id)

if [[ $? != 0 ]]; then
    printf "Failed to fetch analysis output. \n"
else
    # analysis_output_data=$(echo $analysis_output | jq -r ".items[].data[]")
    printf "Fetching analysis output folder id... \n"
    analysis_output_folder_id=$(echo $analysis_output | jq -r ".items[].data[].dataId")

    printf "Downloading analysis output folder with id '$analysis_output_folder_id' to '$local_download_path'... \n"
    analysis_output_download_response=$(icav2 projectdata download $analysis_output_folder_id $local_download_path)
fi

