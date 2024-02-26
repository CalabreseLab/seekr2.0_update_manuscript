#!/bin/bash

#SBATCH -t 24:00:00

# run command: run_downsample_merge.sh /path/to/file1.bam,/path/to/file2.bam,etc experiment_name

# make sure these directories exist, if not then make them.
if [ ! -d downsampled_files ]; then
	mkdir downsampled_files
fi

if [ ! -d merged_downsampled_bamfiles ]; then
	mkdir merged_downsampled_bamfiles
fi
	

file_paths=${1}
experiment_name=${2}
# split file_paths into multiple variables by commas
IFS=',' read -ra file_paths_array <<< "$file_paths"
declare -A alignment_sizes
# for each file path in the array:
for i in "${file_paths_array[@]}"
do
    file_counts=$(samtools view -c $i)
    # store file_counts to a variable
    alignment_sizes["$i"]=${file_counts}

done

# find which of the _counts variables in the smallest
# we will downsample the other files to match the number of counts in this file

for key in "${!alignment_sizes[@]}"; do
  echo "$key: ${alignment_sizes[$key]}"
  # You can perform comparisons here. For instance, to find the max number:
  if [[ -z $min ]] || (( ${alignment_sizes[$key]} < $min )); then
    min=${alignment_sizes[$key]}
    min_key=$key
  fi
done

echo "The file with the smallest alignment size is $min_key with a number $min"


for i in "${file_paths_array[@]}"
do
    file_name=$(echo $i | rev | cut -d'/' -f1 | rev)
    # remove the Aligned.out.sam from the end of the file_name
    file_name=$(echo $file_name | sed 's/Aligned.out.sam//g')

    if [ $i != $min_key ]; then
        percent_diff=$(echo "scale=4; (${alignment_sizes[$i]} - $min) / ${alignment_sizes[$i]}" | bc)
        # calculate 1 - percent_diff
        fraction=$(echo "scale=4; 1 - $percent_diff" | bc)
        # if fraction == 1 then samtools view -b -s $fraction will not work
        if [ $fraction == 1 ] ; then
            samtools view -b $i > downsampled_files/${file_name}_downsampled.bam
        else
            samtools view -b -s $fraction $i > downsampled_files/${file_name}_downsampled.bam
        fi
    else
        samtools view -b $i > downsampled_files/${file_name}_downsampled.bam
    fi

done

merged_file_name=${experiment_name}.bam
full_file_string=""

for i  in "${file_paths_array[@]}"
do
    file_name=$(echo $i | rev | cut -d'/' -f1 | rev)
    # remove the Aligned.out.sam from the end of the file_name
    file_name=$(echo $file_name | sed 's/Aligned.out.sam//g')
    full_file_string="${full_file_string}downsampled_files/${file_name}_downsampled.bam "
samtools merge -f -o merged_downsampled_bamfiles/${merged_file_name} ${full_file_string}
done
