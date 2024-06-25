#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=2-0
#SBATCH --mem=96g

module load samtools


# sort 
samtools sort ENCFF086HIT_RBM15.bam -o ENCFF086HIT_RBM15_sorted.bam

samtools sort ENCFF739LLZ_RBM15.bam -o ENCFF739LLZ_RBM15_sorted.bam

samtools sort ENCFF570EJV_HNRNPK.bam -o ENCFF570EJV_HNRNPK_sorted.bam

samtools sort ENCFF704FID_HNRNPK.bam -o ENCFF704FID_HNRNPK_sorted.bam

samtools sort ENCFF162SAS_MATR3.bam -o ENCFF162SAS_MATR3_sorted.bam

samtools sort ENCFF014KBZ_MATR3.bam -o ENCFF014KBZ_MATR3_sorted.bam

samtools sort ENCFF765BPN_PTBP1.bam -o ENCFF765BPN_PTBP1_sorted.bam

samtools sort ENCFF659RKW_PTBP1.bam -o ENCFF659RKW_PTBP1_sorted.bam

samtools sort ENCFF379LZD_HNRNPM.bam -o EENCFF379LZD_HNRNPM_sorted.bam

samtools sort ENCFF050PTL_HNRNPM.bam -o ENCFF050PTL_HNRNPM_sorted.bam




# index
samtools index ENCFF086HIT_RBM15_sorted.bam

samtools index ENCFF739LLZ_RBM15_sorted.bam

samtools index ENCFF570EJV_HNRNPK_sorted.bam

samtools index ENCFF704FID_HNRNPK_sorted.bam

samtools index ENCFF765BPN_PTBP1_sorted.bam

samtools index ENCFF659RKW_PTBP1_sorted.bam

samtools index EENCFF379LZD_HNRNPM_sorted.bam

samtools index ENCFF050PTL_HNRNPM_sorted.bam

samtools index ENCFF162SAS_MATR3_sorted.bam

samtools index ENCFF014KBZ_MATR3_sorted.bam



# merge replicates of bam files for the same RBP
samtools merge -o RBM15_sorted_merged.bam ENCFF086HIT_RBM15_sorted.bam ENCFF739LLZ_RBM15_sorted.bam

samtools merge -o HNRNPK_sorted_merged.bam ENCFF570EJV_HNRNPK_sorted.bam ENCFF704FID_HNRNPK_sorted.bam

samtools merge -o PTBP1_sorted_merged.bam ENCFF765BPN_PTBP1_sorted.bam ENCFF659RKW_PTBP1_sorted.bam

samtools merge -o HNRNPM_sorted_merged.bam EENCFF379LZD_HNRNPM_sorted.bam ENCFF050PTL_HNRNPM_sorted.bam

samtools merge -o MATR3_sorted_merged.bam ENCFF162SAS_MATR3_sorted.bam ENCFF014KBZ_MATR3_sorted.bam



# index sorted merged bam

samtools index RBM15_sorted_merged.bam

samtools index HNRNPK_sorted_merged.bam

samtools index PTBP1_sorted_merged.bam

samtools index HNRNPM_sorted_merged.bam

samtools index MATR3_sorted_merged.bam



# only keep the second in pair for the paired-ending sequencing results of RBM15, PTBP1, HNRNPM, MATR3

samtools view -h -f 131 PTBP1_sorted_merged.bam > PTBP1_sorted_merged_second.bam

samtools sort PTBP1_sorted_merged_second.bam -o PTBP1_merged_second_sorted.bam

samtools index PTBP1_merged_second_sorted.bam


samtools view -h -f 131 HNRNPM_sorted_merged.bam > HNRNPM_sorted_merged_second.bam

samtools sort HNRNPM_sorted_merged_second.bam -o HNRNPM_merged_second_sorted.bam

samtools index HNRNPM_merged_second_sorted.bam


samtools view -h -f 131 MATR3_sorted_merged.bam > MATR3_sorted_merged_second.bam

samtools sort MATR3_sorted_merged_second.bam -o MATR3_merged_second_sorted.bam

samtools index MATR3_merged_second_sorted.bam


samtools view -h -f 131 RBM15_sorted_merged.bam > RBM15_sorted_merged_second.bam

samtools sort RBM15_sorted_merged_second.bam -o RBM15_merged_second_sorted.bam

samtools index RBM15_merged_second_sorted.bam


