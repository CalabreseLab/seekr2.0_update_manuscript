#!/bin/bash

module load bedtools

sbatch -p general --time=5:00:00 -n 8 -N 1 --mem=96g -o "chunk_rA_multicov.out" --wrap="bedtools multicov -s -D -bams RBM15_merged_second_sorted.bam HNRNPK_sorted_merged.bam MATR3_merged_second_sorted.bam PTBP1_merged_second_sorted.bam HNRNPM_merged_second_sorted.bam -bed ../../xist_sigchunk_rA.bed"

sbatch -p general --time=5:00:00 -n 8 -N 1 --mem=96g -o "control_rA_multicov.out" --wrap="bedtools multicov -s -D -bams RBM15_merged_second_sorted.bam HNRNPK_sorted_merged.bam MATR3_merged_second_sorted.bam PTBP1_merged_second_sorted.bam HNRNPM_merged_second_sorted.bam -bed ../../xist_sigchunk_rA_ck.bed"

sbatch -p general --time=5:00:00 -n 8 -N 1 --mem=96g -o "chunk_rF_multicov.out" --wrap="bedtools multicov -s -D -bams RBM15_merged_second_sorted.bam HNRNPK_sorted_merged.bam MATR3_merged_second_sorted.bam PTBP1_merged_second_sorted.bam HNRNPM_merged_second_sorted.bam -bed ../../xist_sigchunk_rF.bed"

sbatch -p general --time=5:00:00 -n 8 -N 1 --mem=96g -o "control_rF_multicov.out" --wrap="bedtools multicov -s -D -bams RBM15_merged_second_sorted.bam HNRNPK_sorted_merged.bam MATR3_merged_second_sorted.bam PTBP1_merged_second_sorted.bam HNRNPM_merged_second_sorted.bam -bed ../../xist_sigchunk_rF_ck.bed"

sbatch -p general --time=5:00:00 -n 8 -N 1 --mem=96g -o "chunk_rB1_multicov.out" --wrap="bedtools multicov -s -D -bams RBM15_merged_second_sorted.bam HNRNPK_sorted_merged.bam MATR3_merged_second_sorted.bam PTBP1_merged_second_sorted.bam HNRNPM_merged_second_sorted.bam -bed ../../xist_sigchunk_rB1.bed"

sbatch -p general --time=5:00:00 -n 8 -N 1 --mem=96g -o "control_rB1_multicov.out" --wrap="bedtools multicov -s -D -bams RBM15_merged_second_sorted.bam HNRNPK_sorted_merged.bam MATR3_merged_second_sorted.bam PTBP1_merged_second_sorted.bam HNRNPM_merged_second_sorted.bam -bed ../../xist_sigchunk_rB1_ck.bed"

sbatch -p general --time=5:00:00 -n 8 -N 1 --mem=96g -o "chunk_rB2_multicov.out" --wrap="bedtools multicov -s -D -bams RBM15_merged_second_sorted.bam HNRNPK_sorted_merged.bam MATR3_merged_second_sorted.bam PTBP1_merged_second_sorted.bam HNRNPM_merged_second_sorted.bam -bed ../../xist_sigchunk_rB2.bed"

sbatch -p general --time=5:00:00 -n 8 -N 1 --mem=96g -o "control_rB2_multicov.out" --wrap="bedtools multicov -s -D -bams RBM15_merged_second_sorted.bam HNRNPK_sorted_merged.bam MATR3_merged_second_sorted.bam PTBP1_merged_second_sorted.bam HNRNPM_merged_second_sorted.bam -bed ../../xist_sigchunk_rB2_ck.bed"

sbatch -p general --time=5:00:00 -n 8 -N 1 --mem=96g -o "chunk_rD_multicov.out" --wrap="bedtools multicov -s -D -bams RBM15_merged_second_sorted.bam HNRNPK_sorted_merged.bam MATR3_merged_second_sorted.bam PTBP1_merged_second_sorted.bam HNRNPM_merged_second_sorted.bam -bed ../../xist_sigchunk_rD.bed"

sbatch -p general --time=5:00:00 -n 8 -N 1 --mem=96g -o "control_rD_multicov.out" --wrap="bedtools multicov -s -D -bams RBM15_merged_second_sorted.bam HNRNPK_sorted_merged.bam MATR3_merged_second_sorted.bam PTBP1_merged_second_sorted.bam HNRNPM_merged_second_sorted.bam -bed ../../xist_sigchunk_rD_ck.bed"

sbatch -p general --time=5:00:00 -n 8 -N 1 --mem=96g -o "chunk_rE_multicov.out" --wrap="bedtools multicov -s -D -bams RBM15_merged_second_sorted.bam HNRNPK_sorted_merged.bam MATR3_merged_second_sorted.bam PTBP1_merged_second_sorted.bam HNRNPM_merged_second_sorted.bam -bed ../../xist_sigchunk_rE.bed"

sbatch -p general --time=5:00:00 -n 8 -N 1 --mem=96g -o "control_rE_multicov.out" --wrap="bedtools multicov -s -D -bams RBM15_merged_second_sorted.bam HNRNPK_sorted_merged.bam MATR3_merged_second_sorted.bam PTBP1_merged_second_sorted.bam HNRNPM_merged_second_sorted.bam -bed ../../xist_sigchunk_rE_ck.bed"


