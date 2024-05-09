#!/bin/bash

project_name="SGDP"
time_stamp=$(date +"%Y-%m-%d_%H-%M-%S")

printf "Fetching list of projects... \n"

projects=$(icav2 projects list)

printf "$time_stamp: "
if [[ $? != 0 ]]; then
    printf "Failed to fetch list of projects. \n"
else
    project_id=$(echo $projects | jq -r ".items[] | select(.name == \"$project_name\").id")
    if [[ $project_id == "" ]]; then
        printf "No project with name '$project_name' found. \n"
    else
        printf "Found project '$project_name' with ID '$project_id' \n"
    fi
fi
