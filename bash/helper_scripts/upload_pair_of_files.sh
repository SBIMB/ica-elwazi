#!/bin/bash

file_path="$HOME/Documents/sbimb/ica-v2-poc/public/assets/ica_data_uploads/fastq/1_control_rbcLa_2019_minq7"
project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
time_stamp=$(date +"%Y-%m-%d %H:%M:%S")

output_file_prefix=$(echo "${file_path##*/}")
read_1_file=""
read_2_file=""
read_1_analysis_code="read1"
read_2_analysis_code="read2"
ica_path=$(echo "/${file_path##*/}/")

manifest_file="manifest.tsv"
touch $manifest_file

for file in "$file_path"/*; do
    if [[ "$file" == *"_R1_"* ]]; then
        read_1_file=$file
        printf "Read 1 file is '$read_1_file'.\n"
    elif [[ "$file" == *"_R2_"* ]]; then
        read_2_file=$file
        printf "Read 2 file is '$read_2_file'.\n"
    else
        printf "No read files present in directory.\n"
    fi
done

read_1_file_response="read_1_file_response.txt"
read_2_file_response="read_2_file_response.txt"

touch $read_1_file_response
touch $read_2_file_response

printf "[$time_stamp]: "
printf "Uploading read 1 file '$read_1_file'... \n"
read_1_upload_response=$(icav2 projectdata upload $read_1_file $ica_path --project-id $project_id)
echo "$read_1_upload_response" > $read_1_file_response

printf "[$time_stamp]: "
printf "Uploading read 2 file '$read_2_file'... \n"
read_2_upload_response=$(icav2 projectdata upload $read_2_file $ica_path --project-id $project_id)
echo "$read_2_upload_response" > $read_2_file_response

# id of file starts with 'fil.'
read_1_file_id=$(grep -i "\"id\": \"fil\." $read_1_file_response | grep -o 'fil[^\"]*')
read_2_file_id=$(grep -i "\"id\": \"fil\." $read_2_file_response | grep -o 'fil[^\"]*')

printf "[$time_stamp]: "
printf "Writing read file ids to manifest file...\n"

manifest_header_1="file_name"
manifest_header_2="file_id"
manifest_header_3="file_ref"
manifest_tab=$(printf "\t")
read_1_file_ref="$read_1_analysis_code:$read_1_file_id"
read_2_file_ref="$read_2_analysis_code:$read_2_file_id"

cat > $manifest_file << EOF
$manifest_header_1 $manifest_tab $manifest_header_2 $manifest_tab $manifest_header_3
$(echo "${read_1_file##*/}") $manifest_tab $read_1_file_id $manifest_tab $read_1_file_ref
$(echo "${read_2_file##*/}") $manifest_tab $read_2_file_id $manifest_tab $read_2_file_ref
EOF

printf "[$time_stamp]: "
printf "Removing temporary files...\n"
rm $read_1_file_response
rm $read_2_file_response
