#! /bin/bash

#
# Create some empty gvcf.gz files for testing the workflows/prepare_samples pipeline 
#

total_samples=1000
group_size=150

out_dir=test_data

mkdir -p $out_dir

batch_no=1
counter=0

mkdir $out_dir/batch_$batch_no

for s in $(seq 1 $total_samples); do
   echo "sequence_data" >> "$out_dir/batch_$batch_no/sample_$s.gvcf.gz"
   echo "sequence_data index" >> "$out_dir/batch_$batch_no/sample_$s.gvcf.gz.tbi"
   gender=XY
   if [[ $(( RANDOM % 2 )) -eq 0 ]]; then
     gender=XX
   fi
   echo "PLIODY ESTIMATION,,Ploidy Estimation,$gender" > "$out_dir/batch_$batch_no/sample_$s.ploidy_estimation_metrics.csv"
   let counter++
   if [[ $counter -ge $group_size ]]; then
      let batch_no++
      mkdir $out_dir/batch_$batch_no
      counter=0
   fi
done

