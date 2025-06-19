#! /bin/bash
module load sequence

fastq_file_path="/dataA/1000G/SRR1295540_1.fastq.gz"
fastq_file_stats_response_file="fastq_file_stats_response.txt"

fastq_file_stats_response=$(seqkit stats $fastq_file_path)

touch $fastq_file_stats_response_file
echo "$fastq_file_stats_response" > $fastq_file_stats_response_file

fastq_file_sum_len=$(awk '{if ($5); print $5}' $fastq_file_stats_response_file)

echo "$fastq_file_sum_len" > fastq_file_sum_len.txt

fastq_file_size=$(cat fastq_file_sum_len.txt | sed -n '2 p')