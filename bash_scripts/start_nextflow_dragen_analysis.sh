#!/bin/bash

project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
pipeline_id="858708e9-6cdd-42ae-ad97-8b86ee07cdc5"
user_reference="dragen_germline_whole_genome"
storage_size="Medium"
read_1_analysis_code="read1"
read_2_analysis_code="read2"
fastqs_analysis_data_code="fastqs"
fastq_list_data_code="fastq_list",
ref_file_analysis_code="ref_tar"

sample_id="NA12878"
ref_file_id="fil.2ed2474579454e9ff42d08dc4483ed72"
read_1_file_id="fil.19d66b08f8c149a9d3fd08dd9384d2de"
read_2_file_id="fil.b1c8cba850c44b579b6408dd8c9b02f1"
fastq_list_file_id="fil.d27e9cceb0ab495ea60008dd9323ff97"

analysis_response=$(icav2 projectpipelines start nextflow $pipeline_id \
    --user-reference $user_reference \
    --project-id $project_id \
    --storage-size $storage_size \
    --input $ref_file_analysis_code:"$ref_file_id" \
    --input $fastqs_analysis_data_code:"$read_1_file_id,$read_2_file_id" \
    --input $fastq_list_data_code:$fastq_list_file_id \
    --parameters enable_map_align:true \
    --parameters enable_map_align_output:true \
    --parameters output_format:CRAM \
    --parameters enable_duplicate_marking:true \
    --parameters enable_variant_caller:true \
    --parameters vc_emit_ref_confidence:GVCF \
    --parameters vc_enable_vcf_output:true \
    --parameters enable_cnv:true \
    --parameters enable_sv:true \
    --parameters repeat_genotype_enable:false \
    --parameters enable_hla:false \
    --parameters enable_variant_annotation:true \
    --parameters output_file_prefix:"$sample_id")

analysis_response_file="analysis_response.txt"
touch $analysis_response_file
echo "$analysis_response" > $analysis_response_file