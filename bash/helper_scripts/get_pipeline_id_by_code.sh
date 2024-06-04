#!/bin/bash

project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
pipeline_code="basic_pipeline"

pipelines_response=$(icav2 projectpipelines list --project-id $project_id)

pipeline_id=$(echo $pipelines_response | jq -r ".items[].pipeline | select(.code == \"$pipeline_code\").id")

echo $pipeline_id
