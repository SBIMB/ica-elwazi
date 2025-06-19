#!/bin/bash

folder_path="$HOME/Documents/ica_data_uploads/fasta"
folder_name="Citrobacter_30_2_uid32453"
project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
sample_id="4bef7ead-314c-4651-8bb6-da0656dc8332"
sample_name="Citrobacter_30_2_sample"
sample_description="a sample containing Citrobacter_30_2_uid32453 fasta files"
time_stamp=$(date +"%Y-%m-%d_%H:%M:%S")

upload_folder_path="$folder_path/$folder_name"

printf "[$time_stamp]: "
printf "Fetching list of samples... \n"
project_sample_list_response=$(icav2 projectsamples list --project-id $project_id)
if [[ $? != 0 ]]; then
    printf "[$time_stamp]: "
    printf "Failed to fetch list of project samples. \n"
    exit 1
else
    projectsample_id=$(echo $project_sample_list_response | jq -r ".items[].sample | select(.name == \"$sample_name\").id")
    projectsample_status=$(echo $project_sample_list_response | jq -r ".items[].sample | select(.name == \"$sample_name\").id")
    printf "[$time_stamp]: "
    printf "Found project sample '$sample_name' with id '$projectsample_id'. \n"
fi

printf "[$time_stamp]: "
printf "Linking sample '$sample_name' with folder '$sample_name'... \n"

upload_folder_response=$(icav2 projectdata upload $upload_folder_path \
    --project-id $project_id \
    --existing-sample \
    --sample-id $sample_id \
    --sample-name $sample_name)

if [[ $? != 0 ]]; then
    printf "[$time_stamp]: "
    printf "Failed to upload folder '$folder_name' and have it linked to existing sample '$sample_name'. \n"
    exit 1
else
    projectsample_id=$(echo $upload_folder_response | jq -r ".sample.id")
    printf "[$time_stamp]: "
    printf "Created project sample with id '$projectsample_id'. \n"
fi