#!/bin/bash

project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
time_stamp=$(date +"%Y-%m-%d %H:%M:%S")

sample_id="NA12878"
read_1_file_name="NA12878_1.fastq.gz"
read_2_file_name="NA12878_2.fastq.gz"

csv_file="fastq-list-$sample_id.csv"
touch $csv_file

printf "[$time_stamp]: "
printf "Writing to FASTQ list file... \n"
printf "RGID,RGSM,RGLB,Lane,Read1File,Read2File\n" >> $csv_file
printf "$sample_id,$sample_id,RGLB,1,$read_1_file_name,$read_2_file_name\n" >> $csv_file

printf "[$time_stamp]: "
printf "Uploading CSV list '$csv_file'... \n"
csv_file_upload_response=$(icav2 projectdata upload $csv_file --project-id $project_id)
echo "$csv_file_upload_response" > $csv_file_upload_response

# id of file starts with 'fil.'
csv_file_id=$(grep -i "\"id\": \"fil\." $csv_file_upload_response | grep -o 'fil[^\"]*')
