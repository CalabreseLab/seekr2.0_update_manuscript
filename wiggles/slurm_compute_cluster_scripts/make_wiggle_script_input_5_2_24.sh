#!/bin/bash

# run with make_wiggle_script_input_1_30_24.sh </path/to/bed12/files/> </path/to/chrNameLength.txt> <stranded/unstranded> <log10/nolog> <bin size i.e. 50> </path/to/readcounts/ OR 1>

bedfile_path=$1
chrNameLengthfile=$2
strandedness=$3 #stranded or unstranded
log_norm=$4 #y or n
bin_size=$5
readcount_path=$6 #either a path OR the value 1

echo "#!/bin/bash" > run_wiggle_script.sh

for file in $bedfile_path/*.bed
do	
	file_name=$(basename $file ".bed")
	if [ $readcount_path != '1' ]; then
		# get the contents of the corresponding file from readcounts/
		read_counts=$(cat ${readcount_path}/${file_name}_readcounts.txt)
	else
		read_counts='1'
	fi
	
	if [ $strandedness == 'stranded' ]; then
		# identify strand to assign color
		if [[ $file_name == *+* ]]; then
			echo "sbatch --wrap=\"python3 make_wiggle_tracks_1_11_24.py $file $chrNameLengthfile ${file_name} red ${log_norm} ${bin_size} ${read_counts}\"" >> run_wiggle_script.sh
		elif [[ $file_name == *-* ]]; then
			echo "sbatch --wrap=\"python3 make_wiggle_tracks_1_11_24.py $file $chrNameLengthfile ${file_name} blue ${log_norm} ${bin_size} ${read_counts}\"" >> run_wiggle_script.sh
		else
			echo $file_name
			echo "error: strand not clearly detectable"
			exit
		fi
	else
		echo "sbatch --wrap=\"python3 make_wiggle_tracks_1_11_24.py $file $chrNameLengthfile ${file_name} black ${log_norm} ${bin_size} ${read_counts}\"" >> run_wiggle_script.sh
	fi
done
	
