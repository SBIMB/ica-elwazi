#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process batch_files {
	input:
	path search_path
	val batch_size
	val version_no
	val first_batch

	output:
	path 'local-batch-gvcf/*', emit: local_batches
	path 'version-gender/*', emit: genders
	path 'version-batch/*', emit: version_batches

	script:
	"""
	batch=${first_batch} 
	version=${version_no}
	counter=0

	mkdir local-batch-gvcf
	mkdir version-gender
    mkdir version-batch
	
	echo -e "#FID\\tIID\\tSEX" > version-gender/version-\${version}.genders.tsv
	
	for sequence_file in \$(find -L ${search_path} -name "*.gvcf.gz" | xargs readlink -f ); do
		counter=\$((counter+1))
		if [[ ${batch_size} != all && \${counter} -gt ${batch_size} ]]; then
			counter=0
			batch=\$((batch+1))
		fi
		
		echo \${sequence_file} >> local-batch-gvcf/batch-\${batch}.gvcfs.csv
		sample=\${sequence_file%%.*}
		sample_name=\${sample##*/}
		chroms=\$(tail -n1 \${sample}.ploidy_estimation_metrics.csv | cut -f 4 -d , )
		gender=0
		if [ \${chroms} = XX ]; then
          gender=2
		fi
		if [ \${chroms} = XY ]; then
          gender=1
		fi
		echo -e "\${sample_name}\\t\${sample_name}\\t\${gender}" >> version-gender/version-\${version}.genders.tsv

	done

	# We assume that each version should include all previous batches.
    for i in \$(seq \${batch}); do
	    echo \${i} >> version-batch/version-\${version}.batches.txt
    done
	"""
}
