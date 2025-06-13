#!/bin/bash

project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
local_download_path="/home/regan/ica_data_downloads"

array_of_output_folder_ids=("fol.746fc007cb304ac03a8608dd49010d1c" "fol.9814d7855f41483d0c3808dd4847c565" "fol.763b95dd76bf4a9cf12308dd4847c564" "fol.a563d65e340a4a2a6aaf08dd4b28f73b")

for folder_id in ${array_of_output_folder_ids[@]}; 
do
    temp_file="$folder_id.txt"
    touch $temp_file

    printf "Downloading analysis output folder with id '$folder_id' to '$local_download_path'... \n"
    folder_download_response=$(icav2 projectdata download $folder_id $local_download_path --project-id $project_id)
    
    printf "$folder_download_response\n" > $temp_file
done