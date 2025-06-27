#!/bin/bash

set -eo pipefail

icav2=$(command -v icav2 || true)

if [[ x$icav2 == x ]]; then
    echo "*ERROR*: ICA CLI (icav2) not found."
    echo
    echo "- Please download the latest CLI from the following link"
    echo "  https://help.ica.illumina.com/command-line-interface/cli-releasehistory#releases"
    echo "- After download, add the CLI to your \$PATH"
    echo "- Generate ICA API key"
    echo "  - Go to https://<domain>.login.illumina.com"
    echo "  - Click on your user profile"
    echo "  - Select \"Manage API Keys\" and click \"Generate\""
    echo "- Run \"icav2 config set\" to configure CLI (you will need the ICA API key)"
    echo
    echo "Please contact Illumina ICA team if you need help."
    echo
    exit
fi

project_config=$(ls project-data/secret/project-config.*.json | tail -n1)
if [[ x$project_config == x ]]; then
    project_config=project-config.json
fi

if [[ ! -f $project_config ]]; then
    echo "USAGE: "
    echo "  ./2.upload.sh  [project_config_json]"
    echo
    echo "Please generate project config json file using \"1.config.sh \""
    echo
    echo "*ERROR*: Missing project configure json file (default: project-config.json)"
    echo
    exit
fi

key=$(cat $project_config | jq -r .ica_api_key)
project_id=$(cat $project_config | jq -r .ica_job_project.id)
project_name=$(cat $project_config | jq -r .ica_job_project.name)

cd project-data
for d in $(ls | grep -v secret); do
    echo y | icav2 projectdata upload --project-id $project_id -k "$key" $d /meta/
done
cd ..

echo
echo "Successfully uploaded project data to ICA project \"$project_name\""
echo
