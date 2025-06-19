#!/bin/bash

project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
sample_name="Citrobacter_freundii_sample"
sample_description="a sample containing Citrobacter_freundii_4_7_47CFAA fasta files"
time_stamp=$(date +"%Y-%m-%d_%H:%M:%S")

printf "[$time_stamp]: "
printf "Creating project sample '$sample_name'... \n"

projectsample_create_response=$(icav2 projectsamples create $sample_name \
    --project-id $project_id \
    --description $sample_description \
    --user-tag "regan_citrobacter_freundii_sample_01" \
    --technical-tag "citrobacter_freundii_fasta_sample")

sleep 5
if [[ $? == 0 ]]; then
    projectsample_id=$(echo $projectsample_create_response | jq -r ".sample.id")
    printf "[$time_stamp]: "
    printf "Created project sample with id '$projectsample_id'. \n"
else
    printf "[$time_stamp]: "
    printf "Failed to create project sample '$sample_name'. \n"
    exit 1
fi
