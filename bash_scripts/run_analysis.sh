#!/bin/bash

project_id="049307d6-85dd-4cdc-b88d-a740e4e9e550"
pipeline_id="858708e9-6cdd-42ae-ad97-8b86ee07cdc5"
analysis_id="f046c8e2-4add-4716-acd3-c5674b1fed74"
user_reference="regan_dragen_germline_whole_genome_test"
storage_size="Small"
local_download_path="${HOME}/Documents/ica_data_downloads/"
interval_in_seconds=120
sample_id="ERR1019034"
output_directory="/output/ERR1019034/"
read1_file_id="fil.7268e6810bb24d6b978d08dcaaee06d9"
read2_file_id="fil.f1e95d34d124472f19bd08dcaac709b2"
read1_file_name="ERR1019034_1.fastq.gz"
read2_file_name="ERR1019034_2.fastq.gz"
read1_analysis_code="read1:fil.7268e6810bb24d6b978d08dcaaee06d9"
read2_analysis_code="read2:fil.f1e95d34d124472f19bd08dcaac709b2"
reference_analysis_code="ref_tar:fil.2e3fd8d802ee4963da2208dc484ea8f0"
hash_table_config_file="/scratch/reference/hash_table.cfg"

analysis_response=$(icav2 projectpipelines start nextflow $pipeline_id \
--user-reference $user_reference \
--project-id $project_id \
--storage-size $storage_size \
--input $reference_analysis_code \
--input fastqs:$read1_file_id,$read2_file_id \
--parameters enable_map_align:true \
--parameters enable_map_align_output:true \
--parameters output_format:BAM \
--parameters enable_variant_caller:true \
--parameters vc_emit_ref_confidence:BP_RESOLUTION \
--parameters enable_cnv:false \
--parameters enable_sv:true \
--parameters repeat_genotype_enable:true \
--parameters enable_hla:false \
--parameters enable_variant_annotation:true \
--parameters output_file_prefix:"$sample_id")

echo $analysis_response


