#!/bin/bash

project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
pipeline_code="basic_pipeline"
pipeline_id="bfecca03-6443-45bd-b313-e4f555cd0748"
user_reference="regan_dragen_analysis_01"
storage_size="Small"
file_name="NZ_GG704940.fa"
file_id="fil.09b47f75e6014039ac4908dc6dfe2358"
$read_1_analysis_code="read1"
$read_2_analysis_code"read2"
$ref_file_analysis_code"ref_tar"
output_file_directory="/NZ_GG704940/output_files/"
output_file_prefix="NZ_GG704940"

icav2 projectpipelines start nextflow $pipeline_id \
    --user-reference $user_reference \
    --project-id $project_id \
    --storage-size $storage_size \
    --input "$read_1_analysis_code:$read_1_file_id" \
    --input "$read_2_analysis_code:$read_2_file_id" \
    --input "$ref_file_analysis_code:$ref_file_id" \
    --parameters enable-variant-caller:true \
    --parameters RGID:Illumina_RGID \
    --parameters RGSM:$output_file_prefix \
    --parameters output-directory:$output_file_directory \
    --parameters output-file-prefix:$output_file_prefix 
