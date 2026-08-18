#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process make_subshards {
    errorStrategy 'retry'
    maxRetries 3

    input:
    val semaphor
    path config_file
    val version_no
    val shard_ids

    output:
    path "project-data/version-subshard/*"

    script:
    """ 
    function parse_range() {
     input="\$1"
     output=()

     IFS=',' read -ra parts <<< "\$input"
     for part in "\${parts[@]}"; do
 	 if [[ "\$part" == *"-"* ]]; then
       start=\${part%-*}
       end=\${part#*-}
       for ((i=start; i<=end; i++)); do
 	     output+=("\$i")
       done
 	 else
       output+=("\$part")
 	 fi
     done
     printf "%s\n" "\${output[@]}" | sort -n | uniq
    }

    prefix=\$(cat ${config_file} | jq -r .ica_job_project.name | sed 's|-jobs||g')
    key=\$(cat ${config_file} | jq -r .ica_api_key)
    pid=\$(icav2 projects list -k "\$key" | grep \$prefix-msvcf-version-${version_no} | awk '{print \$1}')

    # download dragen.shard.num-subshards.csv
    mkdir -p tmp
    for i in \$(parse_range ${shard_ids}); do
        icav2 projectdata download -k "\$key" --project-id \$pid /data/global-census/shard-\$i/dragen.shard.num-subshards.csv tmp/shard-\$i.dragen.shard.num-subshards.csv
        icav2 projectdata download -k "\$key" --project-id \$pid /data/global-census/shard-\$i/dragen.subshards.stats.tsv tmp/shard-\$i.dragen.subshards.stats.tsv
    done

    mkdir -p project-data/version-subshard
    rm -rf project-data/version-subshard/version-${version_no}.subshards.csv
    rm -rf project-data/version-subshard/version-${version_no}.subshards.stats.tsv
    for i in \$(parse_range ${shard_ids}); do
        cat tmp/shard-\$i.dragen.shard.num-subshards.csv >> project-data/version-subshard/version-${version_no}.subshards.csv
	    cat tmp/shard-\$i.dragen.subshards.stats.tsv >> project-data/version-subshard/version-${version_no}.subshards.stats.tsv
    done
    rm -r tmp

    # upload to job project meta folder
    pid=\$(icav2 projects list -k "\$key" | grep \$prefix-jobs | awk '{print \$1}')
    icav2 projectdata upload -k "\$key" --project-id \$pid project-data/version-subshard/version-${version_no}.subshards.csv /meta/version-subshard/version-${version_no}.subshards.csv
    icav2 projectdata upload -k "\$key" --project-id \$pid project-data/version-subshard/version-${version_no}.subshards.stats.tsv /meta/version-subshard/version-${version_no}.subshards.stats.tsv
 
 """
}
