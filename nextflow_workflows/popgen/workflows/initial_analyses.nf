#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.config_file = './project-data/secret/project-config.json'
params.input_config_dir = './popgen-dragen-igg-demo/project-data/config'
params.gvcf_version = '3.9.5'
params.batch_size = 2
params.version_no = 1
params.samples_dir = './test_data/version-1/'
params.shards = '88,89'
params.chroms = '20'
params.concat_options = 'reduced-site-vcf,reduced-vcf,pgen,site-vcf,vcf'

include { write_config } from '../modules/configure'
include { batch_files } from '../modules/prepare_samples'
include { project_settings ; upload_sequence_files } from '../modules/upload'
include { upload_files as config_uploader } from '../modules/upload'
include { upload_files as batches_uploader } from '../modules/upload'
include { upload_files as genders_uploader } from '../modules/upload'

include { igg_pipeline } from '../workflows/pipeline.nf'


workflow prepare_project {
    take:
    config_file
    input_config
    gvcf_version
    version_no
    batch_no
    batch_size
    samples_dir

    main:
    project_settings(config_file)
    static_config = write_config(input_config, gvcf_version)

    batch_files(samples_dir, batch_size, version_no, batch_no)
    version_batches = batch_files.out.version_batches
    genders = batch_files.out.genders

    upload_sequence_files(
        batch_files.out.local_batches,
        project_settings.out.key,
        project_settings.out.id,
        project_settings.out.name,
        version_no,
        batch_no,
    )

    config_uploader(project_settings.out.key, project_settings.out.id, static_config, '/meta/config/')
    batches_uploader(project_settings.out.key, project_settings.out.id, version_batches, '/meta/version-batch/')
    genders_uploader(project_settings.out.key, project_settings.out.id, genders, '/meta/version-gender/')

    emit:
    batches = upload_sequence_files.out.batch_files
    sequence_dir = upload_sequence_files.out.upload_dir
    static_config
    version_batches
    genders
}

workflow {

    main:

    def config_file = file(params.config_file)
    def input_config = file(params.input_config_dir)
    def gvcf_version = params.gvcf_version
    def version_no = params.version_no
    def batch_no = 1
    def shards = params.shards
    def chroms = params.chroms
    def concat_options = params.concat_options
    def batch_size = params.batch_size
    def samples_dir = file(params.samples_dir)

    prepare_project(
        config_file,
        input_config,
        gvcf_version,
        version_no,
        batch_no,
        batch_size,
        samples_dir,
    )

    static_config = prepare_project.out.static_config.collect()
    batches = prepare_project.out.batches.collect()
    version_batches = prepare_project.out.version_batches.collect()
    version_genders = prepare_project.out.genders.collect()

    igg_pipeline(
        config_file,
        static_config,
        batches,
        prepare_project.out.sequence_dir,
        version_batches,
        version_genders,
        batch_no,
        version_no,
        shards,
        chroms,
        concat_options,
    )

    publish:
    config_file = config_file
    static_config = static_config
    batch_gvcfs = batches
    version_batches = version_batches
    version_genders = version_genders
    census_data = [batch_no: batch_no, files: igg_pipeline.out.census_data]
    msvcf_data = [version_no: version_no, files: igg_pipeline.out.msvcf_data]
}

output {
    config_file {
        path "version-${params.version_no}/project-data/secret/project-config.json"
        mode 'copy'
    }
    static_config {
        path "version-${params.version_no}/project-data/"
        mode 'copy'
    }
    batch_gvcfs {
        path "version-${params.version_no}/project-data/"
        mode 'copy'
    }
    version_batches {
        path "version-${params.version_no}/project-data/"
        mode 'copy'
    }
    version_genders {
        path "version-${params.version_no}/project-data/"
        mode 'copy'
    }
    census_data {
        path { results ->
            "version-${params.version_no}/cohort-census/batch-${results.batch_no}"
        }
        mode 'copy'
    }
    msvcf_data {
        path { results ->
            "version-${params.version_no}/msvcf/version-${results.version_no}"
        }
        mode 'copy'
    }
}
