# Wiggle Track Generation Pipeline
Date: 2/26/24  
Author: Quinn Eberhard  
Lab: Mauro Calabrese  

## Introduction

Welcome to the Wiggle Track generation pipeline. This tutorial will walk you through how to turn your sequencing data into a wiggle script, and if desired, to convert them to bigWigs. Many steps here specify `sbatch` processes, so this is highly recommended for use on a cluster and has not been tested for local machines.

Versions of software used are more for reference than a version-adherence requirement (ex: I have run this with both samtools/1.18 and 1.19, both worked fine). Python 3 is required throughout this process. Please read through all steps even if you do not do them, as some formatting information is important for later steps.




## 1. Align Seq Data
**NECESSARY MODULES/SOFTWARE: star/2.7.10a**

Begin by aligning each replicate. For this example, I use STAR to build the genome index:

```
#!/bin/bash

#SBATCH -t 24:00:00
#SBATCH --mem=150G
#SBATCH -n 8

#STAR --runThreadN 8 --runMode genomeGenerate --genomeDir genomeDir --genomeFastaFiles <genome fasta file here>
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir genomeDir --genomeFastaFiles GRCh38.primary_assembly.genome.fa
```


And once this completes, run the aligment step.

#### If you have **PAIRED-END** data:

```
sbatch -t 24:00:00 -n 24 --mem=50G --error=std_err_%j.err --wrap="STAR --runThreadN 8 --genomeDir genomeDir/ --outFileNamePrefix <path_to_outfile_directory_and_outfile_name> --readFilesIn <path_to_read_1.fastq> <path_to_read_2.fastq>"
```

A real paired-end example:
```
sbatch -t 24:00:00 -n 24 --mem=50G --error=std_err_%j.err --wrap="STAR --runThreadN 8 --genomeDir genomeDir/ --outFileNamePrefix ./alignments/GSE158901_SRR12762158 --readFilesIn GSE158901/SRR12762158_pass_1.fastq GSE158901/SRR12762158_pass_2.fastq"
```

#### If your data is **UNPAIRED**:
```
sbatch -t 24:00:00 -n 24 --mem=50G --error=std_err_%j.err --wrap="STAR --runThreadN 8 --genomeDir genomeDir/ --outFileNamePrefix <path_to_outfile_directory_and_outfile_name> --readFilesIn <path_to_read_1.fastq>"
```






## 2. Downsample and Merge Replicates
**NECESSARY MODULES/SOFTWARE: samtools/1.18**

You may want to downsample and merge your replicates if you plan to generate one wiggle per condition. Downsampling ensures that your replicates contribution equally to the final sample wiggle once you merge the replicates. You may elect to not downsample and merge your replicates if you are interested in identifying potential batch effect or hope to study variation between replicates via the final wiggle tracks.

### If you elect to downsample and merge
You will need to run the `run_downsample_merge_1_30_24.sh` script using the following input format, run this for each set of replicates you wish to downsample:

```
sbatch run_downsample_merge_1_30_24.sh </path/to/file/replicate_1.bam>,</path/to/file/replicate_2.bam> <descriptive name of experiment>
```

If you have more than two replicates, you can add more paths connected by commas as input. If you are running multiple different samples or conditions that each have multiple replicates, all important experiment information such as assay targets or conditions should be specified in your descriptive experiment name.


Run this input file and you will have one BAM file for every set of replicates in the merged_downsampled_bamfiles/ directory. Note that when you run this, you may get errors about:
```
mkdir: cannot create directory ‘downsampled_files’: File exists
mkdir: cannot create directory ‘merged_downsampled_bamfiles’: File exists
```

This is due to the jobs running in tandem, and therefore the intial directory creation step may be flagged multiple times, but this will not affect the actual performance of the job submission, and these errors can be disregarded.




### If you do not elect to downsample and merge


You may want to filter your data - if so run: 
```
bash samtools_sam_to_bam_filter_1_30_24.sh </path/to/samfiles/>
```
If you do not wish to filter your data, please run: 
```
bash samtools_sam_to_bam_1_30_24.sh </path/to/samfiles/>
```






## 3. Split Strands
**NECESSARY MODULE/SOFTWARE: samtools/1.18**

You must now decide whether you want to split your data by strand. If you used a stranded RNA-seq kit (common post ~2016, but not guaranteed, you can look up your particular sequencing kit methods or contact the manufacturer if you are uncertain), it is generally recommended to split by strand for optimal visualization of the data in a wiggle track. 

### If you do not want to split by strand:
You may skip this step. Just note that wherever your final files are stored from the end of [step 2](https://github.com/CalabreseLab/seekr2.0_update_manuscript/blob/main/wiggles/README_general_wiggle_instructions.md#2-downsample-and-merge-replicates) will be the path you input in [step 4](https://github.com/CalabreseLab/seekr2.0_update_manuscript/blob/main/wiggles/README_general_wiggle_instructions.md#4-make-bed12-files).


### If you do want to split by strand:

**A. You must determine whether your data is paired-end or unpaired.** 
If you performed the alignment step yourself, you should already know this, but please confirm whether your data is paired or unpaired if you are using downloaded alignments from online resources. This can be done by studying the SAM flags, located in column 2 in the body of the SAM file. I recommend using this site (https://broadinstitute.github.io/picard/explain-flags.html) to decode your SAM flag values. If the reads are paired, this will be indicated by the checked boxes relevant to the SAM flag property. Please try multiple flags to confirm consistency throughout your alignment.

**B. You must determine whether your data is reverse-stranded.**  
This can also be done from the SAM flags. Using the same website, identify the combinations of boxes that are checked from the various flags in your SAM alignment file.
#### If you have paired-end data:   
						'first in pair' & 'read reverse strand' = 80
						'second in pair' & 'read reverse strand' = 144
						'first in pair' & 'mate reverse strand' = 96
						'second in pair' & 'mate reverse strand' = 160

#### If you have unpaired data:	
As far as I know, this is impossible to determine retroactively without external documentation or your prior processing of the data. If you still do not know, I recommend assuming that your data is NOT reverse stranded as this is decomplicates how you will interpret the SAM flags in post-processing. (With reverse-stranded data you must flip the SAM flag interpretations which can become confusing). If you make the wiggles and discover that your data *is* reverse-stranded, you can simply change the labels then or come back and re-run this for computational reproducability. 


**C. Using the information you obtained above:**  
#### You may wish to filter your data as you split by strand
If so, run:
```
sbatch samtools_filter_binarize_split_strand_1_30_24.sh <paired/unpaired> <forward/reversed> </path/to/bamfiles/>
```
#### You may not wish to filter your data as you split by strand
If so, run:
```
sbatch samtools_binarize_split_strand_1_30_24.sh <paired/unpaired> <forward/reversed> </path/to/bamfiles/>
```
Depending on your previous [step 3](https://github.com/CalabreseLab/seekr2.0_update_manuscript/blob/main/wiggles/README_general_wiggle_instructions.md#3-split-strands) choice, you may have already filtered your data, and could thus run either script here. 






## 4. Make BED12 Files
**NECESSARY MODULES/SOFTWARE: bedtools/2.29**

We must now convert our BAM files into BED12 files to begin making the wiggles. This can be done with:
```
bash bam_to_bed12_1_30_24.sh </path/to/previous/bamfiles>
```

If you split by strand, this will be your path to the strands/ directory. If you did not, then it will be the path to your output from the end of [step 2](https://github.com/CalabreseLab/seekr2.0_update_manuscript/blob/main/wiggles/README_general_wiggle_instructions.md#2-downsample-and-merge-replicates).






## 5. Standardize Signal by Read Counts  
**NECESSARY MODULES/SOFTWARE: samtools/1.18**

If you would like to compare the signal between wiggle tracks (to compare conditions), it is important to standardize by the number of aligned reads in the dataset. If you do wish to standardize your wiggle track signal, run:
```
bash run_count_reads_1_30_24.sh </path/to/bam/files>
```

If you do not wish to do so, please skip this step. 






## 6. Wiggle Script Input
Next, make the input for the wiggle script. There are several things that must be specified about the data, though some of this can be automated for your convenience:


**1. Path and name of each sample's BED12 file**  
**2. The path to the chrNameLength.txt file.** This is generated by STAR during genome indexing and is stored in the genomeDir/ directory if you made the index from [step 1](https://github.com/CalabreseLab/seekr2.0_update_manuscript/blob/main/wiggles/README_general_wiggle_instructions.md#1-align-seq-data), otherwise see bulleted note below.  
**3. Header for output files to be displayed in genome browser**, often short yet descriptive. Use underscores, I recommend specifying if this a +/- stranded, experiment/control dataset in the name along with any other relevant conditions or assay targets.  
**4. Color of the wiggle track** (I like to use separate colors for +/- stranded wiggles for easy visualization)  
**5. Whether to log10 normalize the data**: 'y' for log10 normalization, 'n' for no log10 normalization  
**6. Bin size for the wiggle track** (50 nucleotides is our standard size)  
**7. Number of reads in the sample's alignment from the result of [step 6](https://github.com/CalabreseLab/seekr2.0_update_manuscript/blob/main/wiggles/README_general_wiggle_instructions.md#6-wiggle-script-input) for the sample.** If you skipped [step 5](https://github.com/CalabreseLab/seekr2.0_update_manuscript/blob/main/wiggles/README_general_wiggle_instructions.md#5-standardize-signal-by-read-counts) and do not wish to standardize, put the value 1 here.  

* If you did not use STAR to align your data and do not have access to a STAR genome index for your organism, you can make a chrNameLength.txt by listing the name of each chromosome in "chr#" format on every line, followed by the chromosome length in nucleotides. These should be tab separated.  

### Automatic Input Generation
To generate the serial job submission instructions for each sample automatically, you may be able to run make_wiggle_script_input_1_30_24.sh. This will still require a bit of knowledge about your processing preferences of 1-7 above.
``
sbatch make_wiggle_script_input_1_30_24.sh </path/to/bed12/files/> </path/to/chrNameLength.txt> <stranded/unstranded> <log norm: y or n> <bin size i.e. 50> </path/to/readcounts/ OR 1>
``

Real data example:
```
sbatch make_wiggle_script_input_1_30_24.sh bedfiles/ /proj/seq/data/STAR_genomes_v277/GRCm38_p6_GENCODE_primary/chrNameLength.txt stranded n 50 readcounts/
```

If 'stranded', then colors will be assigned red for + and blue for -. You can manually change the color assignments in this step, if desired. If the 'stranded' parameter is provided, the script expects to find "_+" and "_-" in the file name to determine strandedness and will flag an error otherwise.


Running the above will create the shell script to be used in [step 7](https://github.com/CalabreseLab/seekr2.0_update_manuscript/blob/main/wiggles/README_general_wiggle_instructions.md#7-run-wiggle-script): `run_wiggle_script.sh`

### Manual Input Generation
If you want to create the input file run instructions manually, use the following format for each wiggle track:
```
sbatch --wrap="python3 make_wiggle_tracks_1_11_24.py <bed_file_path/bedfile> <path/to/chrNameLength.txt> <header> <color in R,G,B> <lognormalize y or n> <bin_size> <number of readcount for sample's alignment OR 1>"
```

Real data example:
```
sbatch --wrap="python3 make_wiggle_tracks_1_11_24.py bedfiles/pabpn1_-_if_reverse_stranded_F16.bed /proj/seq/data/STAR_genomes_v277/GRCm38_p6_GENCODE_primary/chrNameLength.txt pabpn1_-_if_reverse_stranded_F16 blue n 50 16157990"
```

You could submit these individually in terminal or by writing the commands line by line into a shell script to submit all jobs at once.






## 7. Run Wiggle Script
To run the `make_wiggle_tracks_1_11_24.py` script, use the previously generated input file. If you automatically generated this in the previous step, simply run:
```
bash  run_wiggle_script.sh
```

The wiggles will develop in the present working directory. Once they have completed you can move them into a new directory for storage. Wiggle generation should only take ~20 minutes, and is very often shorter than that. Be sure to check outputs for any errors.






## 8. Optional: Normalize Wiggle Signal to Controls
If you have a control sample that you would like to normalize your experimental samples to, you can run the following script. This will generate a new wiggle file for each experimental sample that is normalized to the paired control sample. This is an alternative to viewing the experimental wiggle and the control wiggle in the browser together and may provide a more simplified visualization of the data. Your wiggles should have been scaled by read counts.

For each combination of experiment and control sample, run:
```
sbatch --wrap="python3 control_normalize_wiggles_2_20_24.py <path/to/experimental_wiggle.wig> <path/to/control_wiggle.wig>"
```






## 9. Optional: Make bigWigs
**NECESSARY MODULES/SOFTWARE: ucsctools/320**  

Congratulations on generating your wiggle tracks successfully! If you have more than a few and want to be able to upload a view them all relatively quickly, I recommend converting them to bigWigs and using a TrackHub. To convert to bigWigs, run:
```
bash make_bigwigs_1_30_24.sh </path/to/chrNameLength.txt> </path/to/wiggle/files/>
```

Examples:
```
bash make_bigwigs_1_30_24.sh genomeDir/chrNameLength.txt ./wiggle_outputs/
bash make_bigwigs_1_30_24.sh genomeDir/chrNameLength.txt ./
```
The last example would be used if you did not relocate your wiggles into a subdirectory.