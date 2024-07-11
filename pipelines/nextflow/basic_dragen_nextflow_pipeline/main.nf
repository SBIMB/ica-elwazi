nextflow.enable.dsl = 2

process DRAGEN {

    // The container must be a DRAGEN image that is included in an accepted bundle and will determine the DRAGEN version
    container '079623148045.dkr.ecr.us-east-1.amazonaws.com/cp-prod/7ecddc68-f08b-4b43-99b6-aee3cbb34524:latest'
    pod annotation: 'scheduler.illumina.com/presetSize', value: 'fpga-medium'
    pod annotation: 'volumes.illumina.com/scratchSize', value: '1TiB'

    // ICA will upload everything in the "out" folder to cloud storage 
    publishDir 'out', mode: 'symlink'

    input:
        tuple path(read1), path(read2)
        val sample_id
        path ref_tar

    output:
        stdout emit: result
        path '*', emit: output

    script:
        """
        set -ex
        mkdir -p /scratch/reference
        tar -C /scratch/reference -xf ${ref_tar}
        
        /opt/edico/bin/dragen --partial-reconfig HMM --ignore-version-check true
        /opt/edico/bin/dragen --lic-instance-id-location /opt/instance-identity \\
            --output-directory ./ \\
            -1 ${read1} \\
            -2 ${read2} \\
            --intermediate-results-dir /scratch \\
            --output-file-prefix ${sample_id} \\
            --RGID ${sample_id} \\
            --RGSM ${sample_id} \\
            --ref-dir /scratch/reference \\
            --enable-variant-caller true
        """
}

workflow {
    DRAGEN(
        Channel.of([file(params.read1), file(params.read2)]),
        Channel.of(params.sample_id),
        Channel.fromPath(params.ref_tar)
    )
}