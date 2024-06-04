#!/bin/bash

file_name="NZ_GG704942.fa"
file_path="Documents/ica_data_uploads/fasta/NZ_GG704942.fa"
project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
time_stamp=$(date +"%Y-%m-%d_%H:%M:%S")

printf "Checking if file exists at path '$file_path'... \n"

if ! [ -f $file_path ]; then
  printf "File does NOT exist at path '$file_path'. \n"
  exit 1
fi

printf "Uploading file from path '$file_path'... \n"
uploaded_file_response=$(icav2 projectdata upload $file_path --project-id $project_id)

printf "$time_stamp: "

if [[ $? != 0 ]]; then
    printf "Failed to upload file. \n"
else
    printf "$uploaded_file_response\n" > uploaded_file_response.txt
    # id of file starts with 'fil.'
    file_id=$(grep -i "\"id\": \"fil\." uploaded_file_response.txt | grep -o 'fil[^\"]*')
    printf "Uploaded file with ID '$file_id' \n"
fi

rm uploaded_file_response.txt
