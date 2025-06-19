#!/bin/bash

project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
sample_id="NA12878"
time_stamp=$(date +"%Y-%m-%d %H:%M:%S")

output_file_prefix=$sample_id
read_1_file="/spaces/ica/MGI/NA12878_1.fastq.gz"
read_2_file="/spaces/ica/MGI/NA12878_2.fastq.gz"
read_1_analysis_code="read1"
read_2_analysis_code="read2"

read_1_file_response="read_1_file_response.txt"
read_2_file_response="read_2_file_response.txt"

touch $read_1_file_response
touch $read_2_file_response

printf "[$time_stamp]: "
printf "Uploading read 1 file '$read_1_file'... \n"
read_1_upload_response=$(icav2 projectdata upload $read_1_file --project-id $project_id)
echo "$read_1_upload_response" > $read_1_file_response

printf "[$time_stamp]: "
printf "Uploading read 2 file '$read_2_file'... \n"
read_2_upload_response=$(icav2 projectdata upload $read_2_file --project-id $project_id)
echo "$read_2_upload_response" > $read_2_file_response

# id of file starts with 'fil.'
read_1_file_id=$(grep -i "\"id\": \"fil\." $read_1_file_response | grep -o 'fil[^\"]*')
read_2_file_id=$(grep -i "\"id\": \"fil\." $read_2_file_response | grep -o 'fil[^\"]*')

printf "[$time_stamp]: "
printf "Writing read file ids to data file...\n"

data_file="data.txt"

if ! [ -f ${data_file} ]; then
    echo "Data file does not exist. Creating one..."
    touch ${data_file}
fi

printf "[${time_stamp}]: "
printf "Writing file data to existing data file...\n"

printf "sampleId:${sampleId}\n" >> ${data_file}
printf "${read1AnalysisDataCode}:${read_1_file_id}\n" >> ${data_file}
printf "${read2AnalysisDataCode}:${read_2_file_id}\n" >> ${data_file}
printf "read1Name:${read_1_uploaded_file_name}\n" >> ${data_file}
printf "read2Name:${read_2_uploaded_file_name}\n" >> ${data_file}

printf "[$time_stamp]: "
printf "Removing temporary files...\n"
rm $read_1_file_response
rm $read_2_file_response
