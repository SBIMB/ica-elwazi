#!/bin/bash

file_path="$HOME/Documents/ica_data_uploads/fasta/Citrobacter_youngae_ATCC_29220_uid28665/NZ_ABWL00000000.scaffold.fa"
file_name="NZ_GG730303.fa"
project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
pipeline_id="bfecca03-6443-45bd-b313-e4f555cd0748"
analysis_id="2e065006-14c0-4474-b08e-802d1907e75b"
user_reference="regan_test_analysis_03"
storage_size="Small"
local_download_path="${HOME}/Documents/ica_data_downloads/"
interval_in_seconds=120

upload_file_path="$file_path/$file_name"

# We need the fileId and analysisId to start Nextflow analysis
file_id=""

time_stamp=$(date +"%Y-%m-%d %H:%M:%S")
printf "[$time_stamp]: "
printf "Uploading '$file_name'... \n"

upload_file_response=$(icav2 projectdata upload $upload_file_path --project-id $project_id)

if [[ $? != 0 ]]; then
    printf "Failed to upload file '$file_name'. \n"
    exit 1
else
    temp_file_name="$file_name.txt"
    touch $temp_file_name
    printf "$upload_file_response\n" > $temp_file_name
    # id of file starts with 'fil.'
    file_id=$(grep -i "\"id\": \"fil\." $temp_file_name | grep -o 'fil[^\"]*')

    time_stamp=$(date +"%Y-%m-%d %H:%M:%S")
    printf "[$time_stamp]: "
    printf "Uploaded file with ID '$file_id' \n"

    printf "Deleting file '$temp_file_name' \n"
    rm $temp_file_name
fi

printf "Fetching list of project analyses... \n"
projectanalyses_list_response=$(icav2 projectanalyses list --project-id $project_id)

if [[ $? != 0 ]]; then
    printf "Failed to fetch projectanalyses list. \n"
    exit 1
else
    # analysis_id=$(echo $projectanalyses_list_response | jq -r ".items[] | select(.userReference == \"$user_reference\").id")
    printf "Fetching project analysis with id '$analysis_id'... \n"
    analysis_response=$(icav2 projectanalyses input $analysis_id --project-id $project_id)
    if [[ $? != 0 ]]; then
        printf "Failed to fetch projectanalysis with id '$analysis_id'. \n"
        exit 1
    else
        analysis_data_code=$(echo $analysis_response | jq -r ".items[].code")

        printf "Constructing 'analysisCode:fileId' from analysis response and file details... \n"
        file_ref="$analysis_data_code:$file_id"

        printf "Starting Nextflow analysis for input '$file_ref'... \n"
        project_pipeline_start_response=$(icav2 projectpipelines start nextflow $pipeline_id \
            --user-reference $user_reference \
            --project-id $project_id \
            --storage-size $storage_size \
            --input $file_ref)

        printf "Getting project pipeline reference... \n"
        project_pipeline_reference=$(echo $project_pipeline_start_response | jq -r ".reference")  

        printf "Project pipeline reference is '$project_pipeline_reference'.\n"      

        while true;
        do
            printf "Fetching list of project analyses... \n"
            projectanalyses_response=$(icav2 projectanalyses list --project-id $project_id)
            if [[ $? == 0 ]]; then
                printf "Fetching status of project analysis with ID '$analysis_id' and pipeline reference '$project_pipeline_reference'... \n"
                analysis_status=$(echo $projectanalyses_response | jq -r ".items[] | select(.reference == \"$project_pipeline_reference\").status")
                time_stamp=$(date +"%Y-%m-%d %H:%M:%S")
                printf "[$time_stamp]: "
                if [[ $analysis_status == "SUCCEEDED" ]]; then
                    printf "Analysis SUCCEEDED \n"

                    time_stamp=$(date +"%Y-%m-%d %H:%M:%S")
                    printf "[$time_stamp]: "
                    printf "Fetching analysis output of analysis with id '$analysis_id'... \n"
                    analysis_output=$(icav2 projectanalyses output $analysis_id --project-id $project_id)

                    if [[ $? != 0 ]]; then
                        printf "Failed to fetch analysis output. \n"
                    else
                        printf "Fetching analysis output folder id... \n"
                        analysis_output_folder_id=$(echo $analysis_output | jq -r ".items[].data[].dataId")

                        printf "Downloading analysis output folder with id '$analysis_output_folder_id' to '$local_download_path'... \n"
                        analysis_output_download_response=$(icav2 projectdata download $analysis_output_folder_id $local_download_path)
                    fi
                    break;

                elif [[ $analysis_status == "FAILED" ]]; then
                    printf "Analysis FAILED \n"
                    break;

                elif [[ $analysis_status == "FAILED_FINAL" ]]; then
                    printf "Analysis FAILED_FINAL \n"
                    break;

                elif [[ $analysis_status == "ABORTED" ]]; then
                    printf "Analysis ABORTED \n"
                    break;
                else
                    printf "Analysis still in progress... \n"
                fi
            else
                printf "Failed to fetch status of project analysis with pipeline reference '$project_pipeline_reference'... \n"
                break;
            fi
            sleep $interval_in_seconds;
        done
    fi
fi
