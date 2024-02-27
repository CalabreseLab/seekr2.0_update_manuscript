#!/bin/bash

# only binarize the data
# bash samtools_sam_to_bam_1_30_24.sh </path/to/input/sams/> 

if [ ! -d unstranded_unfiltered_bams/ ]; then
	mkdir unstranded_unfiltered_bams
fi

for f in $1/*.sam; do
	if [[ $f == *"Aligned.out.sam"* ]]; then
		# in the event it was aligned with STAR
		new_name=$(basename $f Aligned.out.sam)
	else
		new_name=$(basename $f .sam)
	fi

	samtools view -b $f > unstranded_unfiltered_bams/${new_name}.bam
done

