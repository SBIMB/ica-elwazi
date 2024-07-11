#!/bin/bash

file_path="$HOME/Documents/sbimb/ica-v2-poc/public/assets/ica_data_uploads/fastq/1_control_rbcLa_2019_minq7"
project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
time_stamp=$(date +"%Y-%m-%d_%H:%M:%S")

folder_name=$(echo "${file_path##*/}")
local_upload_path=$(echo "${file_path##*/}/")

mkdir -p $local_upload_path

for file in "$file_path"/*; do
    printf "[$time_stamp]: "
    printf "Copying file '$file' to '$local_upload_path'... \n"
    cp $file $local_upload_path
done

temp_file="upload_folder_response.txt"
touch $temp_file

printf "[$time_stamp]: "
printf "Uploading folder '$local_upload_path'... \n"
upload_folder_response=$(icav2 projectdata upload $local_upload_path --project-id $project_id)
echo "$upload_folder_response" > $temp_file

# id of folder starts with 'fol.'
folder_id=$(grep -i "\"id\": \"fol\." $temp_file | grep -o 'fol[^\"]*')

printf "[$time_stamp]: "
printf "Deleting temporary file...\n"
rm $temp_file

printf "[$time_stamp]: "
printf "Deleting local upload folder...\n"
rm -r $local_upload_path

printf "[$time_stamp]: "
printf "Uploaded folder with ID '$folder_id'. \n"
touch "$folder_name.txt"
echo "$folder_id" > "$folder_name.txt"
