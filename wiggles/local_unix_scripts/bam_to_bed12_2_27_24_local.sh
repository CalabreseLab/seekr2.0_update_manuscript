#!/bin/bash

# run with bash bam_to_bed12_1_30_24.sh </path/to/bam/files>

if [ ! -d bedfiles ]; then
	mkdir bedfiles
fi

for file in $1/*.bam
do
	file_name=$(basename $file ".bam")
	bedtools bamtobed -bed12 -splitD -i ${file} > bedfiles/${file_name}.bed
done
