# run with bash make_wiggle_script_input_10_24_24_local.sh </path/to/bed12/files/> </path/to/chrNameLength.txt> <stranded/unstranded> <log10/nolog as: y/n> <bin size i.e. 50> </path/to/readcounts/ OR 1> <OPTIONAL: normalize read counts by strand: y/n>

bedfile_path=$1
chrNameLengthfile=$2
strandedness=$3 #stranded or unstranded
log_norm=$4 # y or n
bin_size=$5
readcount_path=$6 #either a path OR the value 1
# check if there is a 7th argument for strand normalization
if [[ $# < 7 ]]; then
	norm_by_strand='y'
else
	norm_by_strand=$7 #y or n (default is y if this argument is not provided to match early versions)
fi

echo "#!/bin/bash" > run_wiggle_script.sh

for file in $bedfile_path/*.bed
do	
	file_name=$(basename $file ".bed")
	if [ $readcount_path != '1' ]; then
		if [ ! $norm_by_strand == "y" ]; then
			# remove the "_+_*" and "_-_*" from the file names for matching 
			read_count_file_name=$(echo $file_name | sed 's/_+.*//; s/_-.*//')
			echo "read_count_file_name: $read_count_file_name"
		else
			read_count_file_name=$file_name
		fi
		# get the contents of the corresponding file from readcounts/
		read_counts=$(cat ${readcount_path}/${read_count_file_name}_readcounts.txt)
	else
		read_counts='1'
	fi
	
	if [ $strandedness == 'stranded' ]; then
		# identify strand to assign color
		if [[ $file_name == *+* ]]; then
			echo "python3 make_wiggle_tracks_1_11_24.py $file $chrNameLengthfile ${file_name} red ${log_norm} ${bin_size} ${read_counts}" >> run_wiggle_script.sh
		elif [[ $file_name == *-* ]]; then
			echo "python3 make_wiggle_tracks_1_11_24.py $file $chrNameLengthfile ${file_name} blue ${log_norm} ${bin_size} ${read_counts}" >> run_wiggle_script.sh
		else
			echo $file_name
			echo "error: strand not clearly detectable"
			exit
		fi
	else
		echo "python3 make_wiggle_tracks_1_11_24.py $file $chrNameLengthfile ${file_name} black ${log_norm} ${bin_size} ${read_counts}" >> run_wiggle_script.sh
	fi
done
	
