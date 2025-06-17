#!/bin/bash
file_path="/home/regan/ica_data_uploads/test/upload_me.txt"
project_id="d8c5084a-87fc-429b-b0b6-04bdf1ef1739"
time_stamp=$(date +"%Y-%m-%d %H:%M:%S")

file_name=$(echo "${file_path##*/}")
ica_path="/GVCFs/$file_name"

file_id=""
uploaded_file_path=""

printf "[$time_stamp]: "
printf "Uploading '$file_name' to '$ica_path'... \n"

upload_file_response=$(icav2 projectdata upload $file_path $ica_path --project-id $project_id)
if [[ $? != 0 ]]; then
    printf "Failed to upload file '$file_name'. \n"
    exit 1
else
    temp_file="$file_name.txt"

    touch $temp_file

    printf "$upload_file_response\n" > $temp_file
    # id of file starts with 'fil.'
    file_id=$(grep -i "\"id\": \"fil\." $temp_file | grep -o 'fil[^\"]*')

    printf "[$time_stamp]: "
    printf "Uploaded file with ID '$file_id' \n"
fi

uploaded_file_data_response=$(icav2 projectdata get $file_id)
if [[ $? != 0 ]]; then
    printf "Failed to fetch data about file with id '$file_id'. \n"
    exit 1
else
    uploaded_file_path=$(echo $uploaded_file_data_response | jq -r ".details.path")
    printf "Path of uploaded file is '$uploaded_file_path'. \n"
fi
