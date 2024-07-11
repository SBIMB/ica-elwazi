#!/bin/bash
file_path="$HOME/Documents/sbimb/ica-v2-poc/public/assets/ica_data_uploads/fasta/Yokenella_regensburgei_ATCC_43003_uid65133/NZ_AGCL00000000.scaffold.fa/NZ_JH417871.fa"
project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
time_stamp=$(date +"%Y-%m-%d %H:%M:%S")

file_name=$(echo "${file_path##*/}")
ica_path="/fasta/Yokenella_regensburgei_ATCC_43003_uid65133/NZ_AGCL00000000.scaffold.fa/"

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

manifest_file="manifest.tsv"
manifest_tab=$(printf "\t")

if ! [ -f $manifest_file ]; then
    echo "Manifest file does not exist. Creating one..."
    touch $manifest_file

    manifest_header_1="file_name"
    manifest_header_2="file_id"
    manifest_header_3="file_ica_path"

    printf "[$time_stamp]: "
    printf "Writing file data to manifest...\n"

    printf "$manifest_header_1 $manifest_tab $manifest_header_2 $manifest_tab $manifest_header_3\n" >> $manifest_file
    printf "$file_name $manifest_tab $file_id $manifest_tab $uploaded_file_path\n" >> $manifest_file
else
    printf "[$time_stamp]: "
    printf "Writing file data to existing manifest...\n"

    printf "$file_name $manifest_tab $file_id $manifest_tab $uploaded_file_path\n" >> $manifest_file
fi
