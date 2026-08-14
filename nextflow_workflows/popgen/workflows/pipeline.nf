#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { aggregate_gvcf } from '../modules/igg/step1'
include { wait_for_ica as await_1 } from '../modules/await'
include { aggregate_census } from '../modules/igg/step2'
include { wait_for_ica as await_2 } from '../modules/await'
include { make_subshards } from '../modules/igg/count_subshards'
include { generate_msvcf } from '../modules/igg/step3'
include { wait_for_ica as await_3 } from '../modules/await'
include { concat_msvcf } from '../modules/igg/step4'
include { wait_for_ica as await_4 } from '../modules/await'
include { annotate_variant } from '../modules/igg/step5'
include { wait_for_ica as await_5 } from '../modules/await'

include { results_projects ; cohort_census ; multi_sample_vcf } from '../modules/download'

include { delete_inputs as remove_sequence_data } from '../modules/delete_from_ica'
include { delete_results as remove_census_data } from '../modules/delete_from_ica'
include { delete_results as remove_msvcf_data } from '../modules/delete_from_ica'

workflow igg_pipeline {
    take:
    config_file
    static_config
    batches
    sequence_dir
    version_batches
    version_genders
    batch_no
    version_no
    shards
    chroms
    concat_options

    main:

    step1_jobs = aggregate_gvcf(
        config_file,
        static_config,
        batches,
        version_batches,
        version_genders,
        batch_no,
        shards,
    )

    sync1 = await_1(config_file, step1_jobs).collect()

    remove_sequence_data(sync1, config_file, sequence_dir)

    step2_jobs = aggregate_census(
        sync1,
        config_file,
        static_config,
        batches,
        version_batches,
        version_genders,
        version_no,
        shards,
    )

    sync2 = await_2(config_file, step2_jobs).collect()

    subshards = make_subshards(
        sync2,
        config_file,
        version_no,
        shards,
    ).collect()

    step3_jobs = generate_msvcf(
        sync2,
        config_file,
        static_config,
        batches,
        version_batches,
        version_genders,
        subshards,
        version_no,
        shards,
    )

    sync3 = await_3(config_file, step3_jobs).collect()

    step4_jobs = concat_msvcf(
        sync3,
        config_file,
        static_config,
        batches,
        version_batches,
        version_genders,
        subshards,
        version_no,
        chroms,
        concat_options,
    )

    sync4 = await_4(config_file, step4_jobs).collect()

    step5_jobs = annotate_variant(
        sync4,
        config_file,
        static_config,
        batches,
        version_batches,
        version_genders,
        subshards,
        version_no,
        chroms,
    )

    sync5 = await_5(config_file, step5_jobs).collect()

    results_projects(sync5, config_file, version_no)

    cohort_census(results_projects.out.census_project, batch_no)
    multi_sample_vcf(results_projects.out.msvcf_project, version_no)

    remove_census_data(results_projects.out.census_project, cohort_census.out.src)
    remove_msvcf_data(results_projects.out.msvcf_project, multi_sample_vcf.out.src)

    emit:
    census_data = cohort_census.out.dest
    msvcf_data = multi_sample_vcf.out.dest
}
