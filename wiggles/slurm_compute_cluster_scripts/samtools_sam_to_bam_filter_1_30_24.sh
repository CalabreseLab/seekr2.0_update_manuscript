#!/bin/bash

# only binarize and filter the data, -q 30 is the default (excludes multimappers)
# bash samtools_sam_to_bam_filter_1_30_24.sh </path/to/input/sams/> 

if [ ! -d unstranded_bams/ ]; then
	mkdir unstranded_bams
fi

for f in $1/*.sam; do
	if [[ $f == *"Aligned.out.sam"* ]]; then
		# in the event it was aligned with STAR
		new_name=$(basename $f Aligned.out.sam)
	else
		new_name=$(basename $f .sam)
	fi

	sbatch --wrap="samtools view -b -q 30 $f > unstranded_bams/${new_name}.bam"
done

