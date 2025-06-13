#!/bin/bash

folder_name="dragen_germline_whole_genome-dda020fa-17f8-46b6-b876-ba29ab731007"
folder_ica_path="/dragen_germline_whole_genome-dda020fa-17f8-46b6-b876-ba29ab731007/"
local_download_path="/home/regan/ica_data_downloads/dragen_germline_whole_genome-SRR622459-failed"
project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"

time_stamp=$(date +"%Y-%m-%d_%H-%M-%S")

printf "Downloading folder '$folder_name' from ICA storage at '$folder_ica_path'... \n"
folder_download_response=$(icav2 projectdata download $folder_ica_path $local_download_path --project-id $project_id)
