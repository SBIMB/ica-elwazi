#!/bin/bash

project_id="d8c5084a-87fc-429b-b0b6-04bdf1ef1739"
pipeline_id="993fa5a5-1e3d-43cd-a26b-1ffa053e3c1e"
user_reference="dragen_germline_whole_genome_v4-4"
storage_size="Medium"
cram_analysis_code="crams"
cram_index_analysis_code="crais"
ref_file_analysis_code="ref_tar"

sample_id="LP2100148-DNA_A01"
ref_file_id="fil.62d43552f5fb41e0a50408dd7845ec0e"
cram_file_id="fil.8b815f0b4b7547502ffe08dd931c5122"
cram_index_file_id="fil.6dbef3856a32416d2fff08dd931c5122"

analysis_response=$(icav2 projectpipelines start nextflow $pipeline_id \
    --user-reference $user_reference \
    --project-id $project_id \
    --storage-size $storage_size \
    --input $ref_file_analysis_code:"$ref_file_id" \
    --input $cram_analysis_code:"$cram_file_id" \
    --input $cram_index_analysis_code:"$cram_index_file_id" \
    --parameters enable_map_align:false \
    --parameters enable_map_align_output:false \
    --parameters enable_duplicate_marking:true \
    --parameters enable_variant_caller:true \
    --parameters vc_emit_ref_confidence:GVCF \
    --parameters vc_enable_vcf_output:true \
    --parameters enable_cnv:true \
    --parameters enable_sv:true \
    --parameters repeat_genotype_enable:true \
    --parameters enable_hla:true \
    --parameters enable_variant_annotation:false \
    --parameters output_file_prefix:"$sample_id")

analysis_response_file="cram_analysis_response.txt"
touch $analysis_response_file
echo "$analysis_response" > $analysis_response_file