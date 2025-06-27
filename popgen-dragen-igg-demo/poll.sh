#!/bin/bash

project_config=$(ls project-data/secret/*.json | tail -n1)
key=$(cat $project_config | jq -r .ica_api_key)
project_id=$(cat $project_config | jq -r .ica_job_project.id)

icav2 projectanalyses list -k "$key" --project-id $project_id | grep jobs | awk '{print $NF}' | sort | uniq -c


