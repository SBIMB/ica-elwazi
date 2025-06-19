#!/bin/bash
file_path="/Users/regancannell/Documents/RWCannell/ica-v2-poc/public/assets/ica_data_uploads/fastq/1_control_rbcLa_2019_minq7"
project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
time_stamp=$(date +"%Y-%m-%d %H:%M:%S")

sample_id="1_control_rbcLa_2019_minq7"
read_1_file_name="1_control_rbcLa_2019_minq7_R1_001.fastq"
read_2_file_name="1_control_rbcLa_2019_minq7_R2_001.fastq"
ica_path="/fastq/$sample_id/"

read_1_file=""
read_2_file=""

printf "Looking for read files in '$file_path'...\n"
for file in "$file_path"/*; do
    if [[ "$file" == *"_R1_001.fastq"* ]]; then
        read_1_file=$file
        printf "Read 1 file is '$read_1_file'.\n"
    elif [[ "$file" == *"_R2_001.fastq"* ]]; then
        read_2_file=$file
        printf "Read 2 file is '$read_2_file'.\n"
    else
        printf "No read files present in directory.\n"
    fi
done

printf "[$time_stamp]: "
printf "Uploading '$read_1_file' to '$ica_path'... \n"

read_1_file_response=$(icav2 projectdata upload $read_1_file $ica_path --project-id $project_id)
if [[ $? != 0 ]]; then
    printf "Failed to upload read 1 file '$read_1_file'. \n"
    exit 1
else
    temp_file="$read_1_file.txt"

    touch $temp_file

    printf "$read_1_file_response\n" > $temp_file
    # id of file starts with 'fil.'
    read_1_file_id=$(grep -i "\"id\": \"fil\." $temp_file | grep -o 'fil[^\"]*')

    printf "[$time_stamp]: "
    printf "Uploaded read 1 file with ID '$read_1_file_id' \n"
fi

read_1_file_data_response=$(icav2 projectdata get $read_1_file_id)
if [[ $? != 0 ]]; then
    printf "Failed to fetch data about read 1 file with id '$read_1_file_id'. \n"
    exit 1
else
    read_1_file_path=$(echo $read_1_file_data_response | jq -r ".details.path")
    printf "Path of uploaded read 1 file is '$read_1_file_path'. \n"
fi

printf "[$time_stamp]: "
printf "Uploading '$read_2_file' to '$ica_path'... \n"

read_2_file_response=$(icav2 projectdata upload $read_2_file $ica_path --project-id $project_id)
if [[ $? != 0 ]]; then
    printf "Failed to upload read 2 file '$read_2_file'. \n"
    exit 1
else
    temp_file="$read_2_file.txt"

    touch $temp_file

    printf "$read_2_file_response\n" > $temp_file
    # id of file starts with 'fil.'
    read_2_file_id=$(grep -i "\"id\": \"fil\." $temp_file | grep -o 'fil[^\"]*')

    printf "[$time_stamp]: "
    printf "Uploaded read 2 file with ID '$read_2_file_id' \n"
fi

read_2_file_data_response=$(icav2 projectdata get $read_2_file_id)
if [[ $? != 0 ]]; then
    printf "Failed to fetch data about read 2 file with id '$read_2_file_id'. \n"
    exit 1
else
    read_2_file_path=$(echo $read_2_file_data_response | jq -r ".details.path")
    printf "Path of uploaded read 2 file is '$read_2_file_path'. \n"
fi

csv_file="fastq-list-$sample_id.csv"
touch $csv_file

printf "[$time_stamp]: "
printf "Writing to FASTQ list file... \n"
printf "RGID,RGSM,RGLB,Lane,Read1File,Read2File\n" >> $csv_file
printf "$sample_id,$sample_id,RGLB,1,$read_1_file_name,$read_2_file_name\n" >> $csv_file