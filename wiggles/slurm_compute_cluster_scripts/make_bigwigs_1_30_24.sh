#!/bin/bash

# load ucsctools/320
# run with bash make_bigwigs_1_30_24.sh </path/to/chrNameLength.txt> </path/to/wiggle/files/>

if [ ! -d bigWigoutputs ]; then
	mkdir bigWigoutputs
fi

# first make wiggle track, then run this.
chrNameLengthfile=$1

#wigToBigWig in.wig chrom.sizes out.bw
for file in $2/*.wig
do
	file_name=$(basename $file ".wig")
	sbatch --wrap="wigToBigWig $file ${chrNameLengthfile} bigWigoutputs/${file_name}.bw"
done

