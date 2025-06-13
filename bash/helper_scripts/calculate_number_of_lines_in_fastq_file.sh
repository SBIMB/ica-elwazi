#!/bin/bash

sample_sizes_file="sample_sizes.tsv"
sample_sizes_tab=$(printf "\t")

echo "Sample sizes file does not exist. Creating one..."
touch $sample_sizes_file

sample_sizes_header_1="sample_id"
sample_sizes_header_2="number_of_lines"

printf "[$time_stamp]: "
printf "Writing data to file...\n"

printf "$sample_sizes_header_1 $sample_sizes_tab $sample_sizes_header_2 $sample_sizes_tab $sample_sizes_header_3\n" >> $sample_sizes_file

sample_directory="/dataA/1000G"
for sample in "$sample_directory"/*; do
    if [[ $sample == *.fastq.gz ]]
    then
        sample_id=$(basename ${sample})
        number_of_lines=$(zcat $sample | wc -l)
    fi
    printf "$sample_id $sample_sizes_tab $number_of_lines\n" >> $sample_sizes_file
done
