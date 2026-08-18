#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process generate_msvcf {
    errorStrategy 'retry'
    maxRetries 3
    conda './.conda/popgen'

    input:
    val semaphor
    path config_file
    path static_config, name: "project-data/config/*"
    path batches, name: "project-data/batch-gvcf/*"
    path version_batches, name: "project-data/version-batch/*"
    path genders, name: "project-data/version-gender/*"
    path subshards, name: "project-data/version-subshard/*"
    val version_no
    val shard_ids

    output:
    path "project-analyses/*"

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
   
    for shard_id in \$(parse_range ${shard_ids}); do
        for num_subshards in \$(cat project-data/version-subshard/version-${version_no}.subshards.csv | grep ^\$shard_id, | cut -d ',' -f2); do
            if [[ \$num_subshards == 0 ]]; then continue; fi

            popgen-cli dragen-igg submit --step-name generate-msvcf \
                --input-project-data-folder-path project-data \
                --input-project-config-file-path ${config_file} \
                --output-analysis-json-folder-path project-analyses \
                --version-ids ${version_no} \
                --shard-ids \$shard_id \
                --subshard-ids 1-\$num_subshards \
                --analysis-instance-tier economy
       done
    done
    """
}
