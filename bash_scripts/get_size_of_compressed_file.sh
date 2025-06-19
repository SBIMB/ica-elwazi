#!/bin/bash

file_path="/dataA/1000G/SRR1295540_1.fastq.gz"

compressed_file_size=$(gzip -t $file_path)

if [ -z "$compressed_file_size" ]; then
    echo "Compressed file size is: $compressed_file_size"
    exit 0
else
    echo "Could not get compressed file size: exiting with error."
    exit 1
fi
