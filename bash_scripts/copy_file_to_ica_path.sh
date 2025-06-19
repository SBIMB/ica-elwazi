#!/bin/bash
project_id="d8c5084a-87fc-429b-b0b6-04bdf1ef1739"
time_stamp=$(date +"%Y-%m-%d %H:%M:%S")
initial_ica_file_path="/GVCFs/upload_me.txt"
new_folder="/test"

file_name=$(echo "${initial_ica_file_path##*/}")

file_id="fil.50f2e4e5c18f4eeedbed08ddacd02a05"

printf "[$time_stamp]: "
printf "Copying '$file_name' to folder '$new_folder'... \n"

copy_file_response=$(icav2 projectdata copy $file_id --destination-folder $new_folder --source-project-id $project_id --project-id $project_id)

printf "$copy_file_response\n"