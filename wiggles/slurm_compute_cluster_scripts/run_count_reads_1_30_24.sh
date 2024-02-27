#!/bin/bash

# run with bash run_count_reads_1_30_24.sh </path/to/bam/files>
if [ ! -d readcounts ]; then
	mkdir readcounts
fi

for file in $1/*.bam
do
	file_name=$(basename $file ".bam")
	sbatch --wrap="samtools view -c $file > readcounts/${file_name}_readcounts.txt"
done
