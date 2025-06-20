#!/bin/bash

project_id="d8c5084a-87fc-429b-b0b6-04bdf1ef1739"
pipeline_id="05b7e493-89d4-445d-853b-83addad39cca"
user_reference="dragen_joint_pedigree_calling_v4-4-4"
storage_size="Medium"
cram_analysis_code="crams"
cram_index_analysis_code="crais"
ref_file_analysis_code="ref_tar"
tn_tsv_files_analysis_code="tn_tsvs"
gvcf_files_analysis_code="gvcfs"

sample_id="LP2100148-DNA_E01"
ref_file_id="fil.62d43552f5fb41e0a50408dd7845ec0e"
cram_file_id="fil.3d563580671047883ad708ddae271f1f"
cram_index_file_id="fil.280b6de833fa4a1b481908ddae1f44ab"

tn_tsv_file_1_id="fil.2ec1d90df8f14c00039008ddacb1f1e0"
tn_tsv_file_2_id="fil.a13b5d66f2e44188e34908ddae03c1a8"
tn_tsv_file_3_id="fil.060e5f1bdacc42b841fa08ddae1f75f3"
tn_tsv_file_4_id="fil.3c700988e8044a18c47b08ddacb4dc92"

gvcf_file_1_id="fil.6fcbfa8cb72e4061038a08ddacb1f1e0"
gvcf_file_2_id="fil.48898b6b0f3c47773a0008ddae271f1f"
gvcf_file_3_id="fil.7f803bdba2e14c0841f608ddae1f75f3"
gvcf_file_4_id="fil.7f30fdd23cc6474b475d08ddae1f44ab"

analysis_response=$(icav2 projectpipelines start nextflow $pipeline_id \
    --user-reference $user_reference \
    --project-id $project_id \
    --storage-size $storage_size \
    --input $ref_file_analysis_code:"$ref_file_id" \
    --input $cram_analysis_code:"$cram_file_id" \
    --input $cram_index_analysis_code:"$cram_index_file_id" \
    --input $tn_tsv_files_analysis_code:"$tn_tsv_file_1_id,$tn_tsv_file_2_id,$tn_tsv_file_3_id,$tn_tsv_file_4_id" \
    --input $gvcf_files_analysis_code:"$gvcf_file_1_id,$gvcf_file_2_id,$gvcf_file_3_id,$gvcf_file_4_id" \
    --parameters enable_joint_genotyping:true \
    --parameters enable_variant_caller:true \
    --parameters vc_emit_ref_confidence:GVCF \
    --parameters vc_enable_vcf_output:true \
    --parameters enable_cnv:true \
    --parameters enable_sv:true \
    --parameters repeat_genotype_enable:false \
    --parameters enable_map_align:false \
    --parameters enable_map_align_output:false \
    --parameters enable_duplicate_marking:false \
    --parameters enable_hla:false \
    --parameters enable_variant_annotation:false \
    --parameters output_file_prefix:"$sample_id")

analysis_response_file="joint_genotyping_analysis_response.txt"
touch $analysis_response_file
echo "$analysis_response" > $analysis_response_file