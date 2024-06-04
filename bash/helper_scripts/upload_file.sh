#!/bin/bash

file_path="$HOME/Documents/ica_data_uploads/fasta/Yokenella_regensburgei_ATCC_43003_uid65133/NZ_AGCL00000000.scaffold.fa"
file_name="NZ_JH417859.fa"
project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
time_stamp=$(date +"%Y-%m-%d_%H:%M:%S")

upload_file_path="$file_path/$file_name"

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

    printf "[$time_stamp]: "
    printf "Uploaded file with ID '$file_id' \n"
fi

rm $temp_file_name