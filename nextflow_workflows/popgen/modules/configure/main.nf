#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process write_config {
    input:
    path input_config, name: 'inputs'
    val gvcf_version

    output:
    path 'config/*'

    script:
    """
    # Copy existing config files, the data in the demo is suitable
    mkdir config
    for c in annodb.txt chroms.txt low-conf.bed ref.csv shards.txt; do
      cp ${input_config}/\${c} config/
    done

    # Create the gvcf version config file
    echo ${gvcf_version} > config/gvcf-version.txt
  """
}
