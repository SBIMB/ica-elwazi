#!/bin/bash

file_path="/home/regan/ica_data_uploads/AGenDA/LP2100148-DNA_A01"
project_id="d8c5084a-87fc-429b-b0b6-04bdf1ef1739"
pipeline_id="993fa5a5-1e3d-43cd-a26b-1ffa053e3c1e"
reference_file_id="fil.62d43552f5fb41e0a50408dd7845ec0e"
sample_id="V350117642_L02"
time_stamp=$(date +"%Y-%m-%d %H:%M:%S")

cram_file_analysis_data_code="cram"
cram_index_file_analysis_data_code="cram_index"
cram_file="/home/regan/ica_data_uploads/AGenDA/LP2100148-DNA_A01/LP2100148-DNA_A01.cram"
cram_index_file="/home/regan/ica_data_uploads/AGenDA/LP2100148-DNA_A01/LP2100148-DNA_A01.cram.crai"
cram_file_analysis_code="cram"
cram_index_file_analysis_code="cram_index"

cram_file_response="cram_file_response.txt"
cram_index_file_response="cram_index_file_response.txt"

touch $cram_file_response
touch $cram_index_file_response

printf "[$time_stamp]: "
printf "Uploading cram file '$cram_file'... \n"
cram_file_upload_response=$(icav2 projectdata upload $cram_file --project-id $project_id)
echo "$cram_file_upload_response" > $cram_file_response

printf "[$time_stamp]: "
printf "Uploading cram index file '$cram_index_file'... \n"
cram_index_file_upload_response=$(icav2 projectdata upload $cram_index_file --project-id $project_id)
echo "$cram_index_file_upload_response" > $cram_index_file_response

# id of file starts with 'fil.'
cram_file_id=$(grep -i "\"id\": \"fil\." $cram_file_response | grep -o 'fil[^\"]*')
cram_index_file_id=$(grep -i "\"id\": \"fil\." $cram_index_file_response | grep -o 'fil[^\"]*')

printf "[$time_stamp]: "
printf "Writing file ids to data file...\n"

data_file="data.txt"

if ! [ -f ${data_file} ]; then
    echo "Data file does not exist. Creating one..."
    touch ${data_file}
fi

printf "[${time_stamp}]: "
printf "Writing file data to existing data file...\n"

printf "sampleId:$sample_id\n" >> $data_file
printf "$cram_file_analysis_code:$cram_file_id\n" >> $data_file
printf "$cram_index_file_analysis_code:$cram_index_file_id\n" >> $data_file

printf "[$time_stamp]: "
printf "Removing temporary files...\n"
rm $cram_file_response
rm $cram_index_file_response
