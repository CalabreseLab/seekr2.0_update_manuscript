#!/bin/bash

# binarizes to a bam file, and splits by strand
# bash samtools_filter_binarize_split_strands_1_30_24.sh <paired/unpaired> <forward/reversed> </path/to/samfiles/>
if [ ! -d strands ]; then
	mkdir strands
fi


paired_end=$1
#echo $paired_end
reversed_seq=$2
#echo $reversed_seq
samfile_path=$3


if [ $paired_end == "unpaired" ]; then
	

	if [ $reversed_seq == "forward" ]; then

		for file in $samfile_path/*.sam; do

			file_name=$(basename $file .sam)
			# since this is forward-stranded data, SAM flags can be read literally
			# Here we map to the "reverse" strand, forcing the 16 flag
			samtools view -b -f 16 $file > strands/"$file_name"_-_if_forward_stranded_f16.bam

			# Next we map to the "forward" strand, excluding the 16 flag
			samtools view -b -F 16 $file > strands/"$file_name"_+_if_forward_stranded_F16.bam
		done




	elif [ $reversed_seq == "reversed" ]; then

		for file in $samfile_path/*.sam; do

			file_name=$(basename $file .sam)
                        # since this is reverse-stranded data, SAM flags mean the inverse
                        # Here we map to the "reverse" strand, forcing the 16 flag
						# However, for our reversed data, this is the forward strand
                        samtools view -b -f 16 $file > strands/"$file_name"_+_if_reverse_stranded_f16.bam

                        # Next we map to the "forward" strand, excluding the 16 flag
						# Which is actually the reverse strand for our reversed data
                        samtools view -b -F 16 $file > strands/"$file_name"_-_if_reverse_stranded_F16.bam
                done



	else
		echo "Please use the following format to run this script: bash samtools_filter_binarize_split_strands_1_30_24.sh <paired/unpaired> <forward/reversed> </path/to/samfiles.sam>"
		exit

	fi	
	


elif [ $paired_end == "paired" ]; then
	

	if [ $reversed_seq == "forward" ]; then
		

		for file in $samfile_path/*.sam; do
			file_name=$(basename $file .sam)
			# since this is forward-stranded data, SAM flags can be read literally
			# Here we map to the "reverse" strand, forcing the 83 and 163 flags
			samtools view -b -f 83 $file > strands/"$file_name"_-_83.bam
			samtools view -b -f 163 $file > strands/"$file_name"_-_163.bam
			samtools merge -f strands/"$file_name"_-_if_forward_stranded_f83_f163.bam strands/"$file_name"_-_83.bam strands/"$file_name"_-_163.bam
			rm strands/"$file_name"_-_83.bam strands/"$file_name"_-_163.bam

			# Next we map to the "forward strand", forcing the 99 and 147 flags
			samtools view -b -f 99 $file > strands/"$file_name"_+_99.bam
			samtools view -b -f 147 $file > strands/"$file_name"_+_147.bam
			samtools merge -f strands/"$file_name"_+_if_forward_stranded_f99_f147.bam strands/"$file_name"_+_99.bam strands/"$file_name"_+_147.bam
			rm strands/"$file_name"_+_99.bam strands/"$file_name"_+_147.bam
		done




	elif [ $reversed_seq == "reversed" ]; then
		
		for file in $samfile_path/*.sam; do

			file_name=$(basename $file .sam)
			# since this is reverse-stranded data, SAM flags mean the inverse
			# Here we map to the "reverse" strand, forcing the 83 and 163 flags
			# However, for our reversed data, this is the forward strand
			samtools view -b -f 83 $file > strands/"$file_name"_+_83.bam
			samtools view -b -f 163 $file > strands/"$file_name"_+_163.bam
			samtools merge -f strands/"$file_name"_+_if_reverse_stranded_f83_f163.bam strands/"$file_name"_+_83.bam strands/"$file_name"_+_163.bam
			rm strands/"$file_name"_+_83.bam strands/"$file_name"_+_163.bam

			# Next we map to the "forward" strand, forcing the 99 and 147 flag
			# Which is actually the reverse strand for our reversed data
			samtools view -b -f 99 $file > strands/"$file_name"_-_99.bam
			samtools view -b -f 147 $file > strands/"$file_name"_-_147.bam
			samtools merge -f strands/"$file_name"_-_if_reverse_stranded_f99_f147.bam strands/"$file_name"_-_99.bam strands/"$file_name"_-_147.bam
			rm strands/"$file_name"_-_99.bam strands/"$file_name"_-_147.bam
		done



	else
		echo "Please use the following format to run this script: bash samtools_filter_binarize_split_strands_1_30_24.sh <paired/unpaired> <forward/reversed> </path/to/samfiles/>"
		exit
	fi
	
	
	
else
	echo "Please use the following format to run this script: bash samtools_filter_binarize_split_strands_1_30_24.sh <paired/unpaired> <forward/reversed> </path/to/samfiles/>"
	exit
fi
