#!/bin/bash

input_project_data_folder_path=project-data
output_analysis_json_folder_path=project-analyses
input_project_config_file_path=$(ls project-data/secret/project-config.*.json | tail -n1)
analysis_instance_tier=economy

function parse_range() {

    input="$1"
    output=()

    IFS=',' read -ra parts <<< "$input"
    for part in "${parts[@]}"; do
	if [[ "$part" == *"-"* ]]; then
            start=${part%-*}
            end=${part#*-}
            for ((i=start; i<=end; i++)); do
		output+=("$i")
            done
	else
            output+=("$part")
	fi
    done
    printf "%s\n" "${output[@]}" | sort -n | uniq

}

function step1 () {
    # aggregate gvcf files

    batch_ids=$1
    shard_ids=$2
    step_name="aggregate-gvcf"

    popgen-cli dragen-igg submit \
        --input-project-data-folder-path $input_project_data_folder_path \
        --input-project-config-file-path $input_project_config_file_path \
        --output-analysis-json-folder-path $output_analysis_json_folder_path \
        --step-name $step_name \
        --batch-ids $batch_ids \
        --shard-ids $shard_ids \
        --track-variants false \
        --analysis-instance-tier $analysis_instance_tier
}

function step2 () {
    # aggregate census files

    version_id=$1
    shard_ids=$2
    step_name="aggregate-census"

    popgen-cli dragen-igg submit \
        --input-project-data-folder-path $input_project_data_folder_path \
        --input-project-config-file-path $input_project_config_file_path \
        --output-analysis-json-folder-path $output_analysis_json_folder_path \
        --step-name $step_name \
        --version-ids $version_id \
        --shard-ids $shard_ids \
        --analysis-instance-tier $analysis_instance_tier
}

function count_subshards () {
    # get number of subshards per shard and upload to ica job project
    # note:
    # - run this before step3 (so that the number of subshard per shard is updated)
    # - this part of demo will be simplified as PopGen CLI utils

    version_id=$1
    shard_ids=$2

    prefix=$(cat $input_project_config_file_path | jq -r .ica_job_project.name | sed 's|-jobs||g')
    key=$(cat $input_project_config_file_path | jq -r .ica_api_key)
    pid=$(icav2 projects list -k "$key" | grep $prefix-msvcf-version-$version_id | awk '{print $1}')

    # download dragen.shard.num-subshards.csv
    mkdir -p tmp
    for i in $(parse_range $shard_ids); do
        icav2 projectdata download -k "$key" --project-id $pid /data/global-census/shard-$i/dragen.shard.num-subshards.csv tmp/shard-$i.dragen.shard.num-subshards.csv
        icav2 projectdata download -k "$key" --project-id $pid /data/global-census/shard-$i/dragen.subshards.stats.tsv tmp/shard-$i.dragen.subshards.stats.tsv
    done

    mkdir -p project-data/version-subshard
    rm -rf project-data/version-subshard/version-$version_id.subshards.csv
    rm -rf project-data/version-subshard/version-$version_id.subshards.stats.tsv
    for i in $(parse_range $shard_ids); do
        cat tmp/shard-$i.dragen.shard.num-subshards.csv >> project-data/version-subshard/version-$version_id.subshards.csv
	    cat tmp/shard-$i.dragen.subshards.stats.tsv >> project-data/version-subshard/version-$version_id.subshards.stats.tsv
    done
    rm -r tmp

    echo "shard #subshard list written to project-data/version-subshard/version-$version_id.subshards.csv"
    echo "shard #subshard stats written to project-data/version-subshard/version-$version_id.subshards.stats.tsv"

    # upload to job project meta folder
    pid=$(icav2 projects list -k "$key" | grep $prefix-jobs | awk '{print $1}')
    icav2 projectdata upload -k "$key" --project-id $pid project-data/version-subshard/version-$version_id.subshards.csv /meta/version-subshard/version-$version_id.subshards.csv
    icav2 projectdata upload -k "$key" --project-id $pid project-data/version-subshard/version-$version_id.subshards.stats.tsv /meta/version-subshard/version-$version_id.subshards.stats.tsv
}

function step3 () {
    # generat msVCF, run ML filtering, generate reduced msVCF, and PGEN files
    version_id=$1
    shard_ids=$2
    step_name="generate-msvcf"

    for shard_id in $(parse_range $shard_ids); do
        for num_subshards in $(cat project-data/version-subshard/version-$version_id.subshards.csv | grep ^$shard_id, | cut -d\, -f2); do
            if [[ $num_subshards == 0 ]]; then continue; fi
            popgen-cli dragen-igg submit \
                --input-project-data-folder-path $input_project_data_folder_path \
                --input-project-config-file-path $input_project_config_file_path \
                --output-analysis-json-folder-path $output_analysis_json_folder_path \
                --step-name $step_name \
                --version-ids $version_id \
                --shard-ids $shard_id \
                --subshard-ids 1-$num_subshards \
                --analysis-instance-tier $analysis_instance_tier
        done
    done
}

function step4 () {
    # concat per output file format
    # for demo, 5 file formats are grouped into 1 job

    version_id=$1
    chrom_ids=$2
    step_name="concat-msvcf"

    for concat_options in reduced-site-vcf,reduced-vcf,pgen,site-vcf,vcf; do
        popgen-cli dragen-igg submit \
            --input-project-data-folder-path $input_project_data_folder_path \
            --input-project-config-file-path $input_project_config_file_path \
            --output-analysis-json-folder-path $output_analysis_json_folder_path \
            --step-name $step_name \
            --version-ids $version_id \
            --chrom-ids $chrom_ids \
            --concat-options $concat_options \
            --analysis-instance-tier $analysis_instance_tier
    done
}

function step5 () {
    # annotate variant

    version_id=$1
    chrom_ids=$2
    step_name="annotate-variant"

    popgen-cli dragen-igg submit \
        --input-project-data-folder-path $input_project_data_folder_path \
        --input-project-config-file-path $input_project_config_file_path \
        --output-analysis-json-folder-path $output_analysis_json_folder_path \
        --step-name $step_name \
        --version-ids $version_id \
        --chrom-ids $chrom_ids \
        --analysis-instance-tier $analysis_instance_tier
}

$@

# ------- demo version1 (batch1) -------

# version_id=1
# batch_ids=1
# shard_ids=88-89
# chrom_ids=20
# step1 $batch_ids $shard_ids
# step2 $version_id $shard_ids
# count_subshards $version_id $shard_ids
# step3 $version_id $shard_ids
# step4 $version_id $chrom_ids
# step5 $version_id $chrom_ids

# ------- demo version2 (batch1 + batch2) -------

# version_id=2
# batch_ids=2
# shard_ids=88-89
# chrom_ids=20
# step1 $batch_ids $shard_ids
# step2 $version_id $shard_ids
# count_subshards $version_id $shard_ids
# step3 $version_id $shard_ids
# step4 $version_id $chrom_ids
# step5 $version_id $chrom_ids
