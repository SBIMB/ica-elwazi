#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.config_file = './project-data/secret/project-config.json'
params.previous_results_dir = './results/'
params.version = 2
params.batch = params.version
params.batch_size = 2
params.samples_dir = './test_data/version-2'
params.shards = '88,89'
params.chroms = '20'
params.concat_options = 'reduced-site-vcf,reduced-vcf,pgen,site-vcf,vcf'

include { batch_files } from '../modules/prepare_samples'
include { results_projects } from '../modules/download'
include { project_settings ; upload_sequence_files } from '../modules/upload'
include { upload_directory as census_uploader } from '../modules/upload'
include { upload_files as config_uploader } from '../modules/upload'
include { upload_files as batches_uploader } from '../modules/upload'
include { upload_files as genders_uploader } from '../modules/upload'

include { igg_pipeline } from '../workflows/pipeline.nf'

workflow prepare_project {
    take:
    config_file
    static_config
    version_no
    batch_no
    batch_size
    samples_dir

    main:
    project_settings(config_file)

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

workflow restore_outputs {
    take:
    config_file
    results_dir

    main:
    // version_no = 0. cohort-census should use the same project for each version.
    results_projects('ok', config_file, 0)
    project_settings(config_file)

    files = results_dir / '*' / 'cohort-census' / '*'

    census_uploader(
        project_settings.out.key,
        results_projects.out.census_project,
        files,
        '*/cohort-census/*',
    )

    emit:
    census_project = results_projects.out.census_project
}

workflow {

    main:

    def config_file = file(params.config_file)
    def previous_data = file(params.previous_results_dir)
    def version_no = params.version
    def batch_no = params.batch
    def shards = params.shards
    def chroms = params.chroms
    def concat_options = params.concat_options
    def batch_size = params.batch_size
    def samples_dir = file(params.samples_dir)

    static_config = channel.fromPath(previous_data / 'project-data/config/*').collect()

    prior_batches = channel.fromPath(previous_data / 'project-data/batch-gvcf/*').collect()

    prepare_project(
        config_file,
        static_config,
        version_no,
        batch_no,
        batch_size,
        samples_dir,
    )

    restore_outputs(config_file, previous_data)

    batches = prepare_project.out.batches.mix(prior_batches).collect()
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
        path 'project-data/secret/project-config.json'
        mode 'copy'
    }
    static_config {
        path 'project-data/'
        mode 'copy'
    }
    batch_gvcfs {
        path 'project-data/'
        mode 'copy'
    }
    version_batches {
        path 'project-data/'
        mode 'copy'
    }
    version_genders {
        path 'project-data/'
        mode 'copy'
    }
    census_data {
        path { results ->
            "cohort-census/${results.batch_no}"
        }
        mode 'copy'
    }
    msvcf_data {
        path { results ->
            "msvcf/${results.version_no}"
        }
        mode 'copy'
    }
}
