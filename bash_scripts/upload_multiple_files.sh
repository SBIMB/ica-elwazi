#!/bin/bash
project_id="d8c5084a-87fc-429b-b0b6-04bdf1ef1739"
time_stamp=$(date +"%Y-%m-%d %H:%M:%S")

# path of text file with list of GVCFs
list_of_files="files_to_be_uploaded.txt"

# Check if the file exists
if [[ ! -f "$list_of_files" ]]; then
  echo "Error: File '$list_of_files' not found."
  exit 1
fi

# Loop through each line in the file
while IFS= read -r file_path; do
    printf "[$time_stamp]: "
    printf "Checking file path: $file_path\n"

    file_name=$(echo "${file_path##*/}")

    upload_file_response=$(icav2 projectdata upload $file_path --project-id $project_id)
    if [[ $? != 0 ]]; then
        printf "[$time_stamp]: "
        printf "Failed to upload file '$file_name'. \n"
        exit 1
    else
        printf "[$time_stamp]: "
        printf "Successfully uploaded file '$file_name' to ICA platform. \n"
    fi
done < "$list_of_files"
