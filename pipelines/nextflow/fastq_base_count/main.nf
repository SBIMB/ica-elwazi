#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process getFastqFileSize {
    debug true
    input:
    tuple val(sampleId), file(reads)

    output:
    path "data.txt", emit: dataFile

    script:
    def (read_1_file, read_2_file) = reads
    """
    #!/bin/bash
    time_stamp=\$(date +"%Y-%m-%d %H:%M:%S")
    read_1_size=""
    read_2_size=""

    read_1_stats_response_file="read_1_stats_response.tsv"
    read_2_stats_response_file="read_2_stats_response.tsv"

    touch \${read_1_stats_response_file}
    touch \${read_2_stats_response_file}

    printf "[\${time_stamp}]: "
    printf "Getting stats of read 1 file '${read_1_file}'... \n"
    read_1_stats_response=\$(seqkit stats ${read_1_file})
    echo "\${read_1_stats_response}" > \${read_1_stats_response_file}

    printf "[\${time_stamp}]: "
    printf "Getting stats of read 2 file '${read_2_file}'... \n"
    read_2_stats_response=\$(seqkit stats ${read_2_file})
    echo "\${read_2_stats_response}" > \${read_2_stats_response_file}

    read_1_sum_len=\$(awk '{if (\$5); print \$5}' \${read_1_stats_response_file})
    read_2_sum_len=\$(awk '{if (\$5); print \$5}' \${read_2_stats_response_file})

    echo "\${read_1_sum_len}" > read_1_sum_len.txt
    echo "\${read_2_sum_len}" > read_2_sum_len.txt

    read_1_size=\$(cat read_1_sum_len.txt | sed -n '2 p')
    read_2_size=\$(cat read_2_sum_len.txt | sed -n '2 p')
    
    data_file="data.txt"

    if ! [ -f \${data_file} ]; then
        echo "Data file does not exist. Creating one..."
        touch \${data_file}
    fi

    printf "[\${time_stamp}]: "
    printf "Writing file data to existing data file...\n"

    printf "sampleId:${sampleId}\n" >> \${data_file}
    printf "read_1_size:\${read_1_size}\n" >> \${data_file}
    printf "read_2_size:\${read_2_size}\n" >> \${data_file}
    """
}

process printData {
    debug true

    input:
    path(dataFile)

    output:
    stdout

    script:
    """
    timeStamp=\$(date +"%Y-%m-%d %H:%M:%S")
    printf "[\${timeStamp}]: Printing contents of dataFile...\n"
    cat ${dataFile}
    """
}

workflow {
    readsPairFilesUploadPath = params.readsPairFilesUploadPath
    fastqFilePairs = Channel.fromFilePairs(readsPairFilesUploadPath, checkIfExists:true)

    getFastqFileSize(fastqFilePairs)
    printData(getFastqFileSize.out.dataFile)
}
