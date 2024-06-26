#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3-0
#SBATCH --mem=160g



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -fm -pb -o v43_full_1

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_full_1.csv -o v43vXIST_full_1



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -fm -pb -o v43_full_2

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_full_2.csv -o v43vXIST_full_2



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -fm -pb -o v43_full_3

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_full_3.csv -o v43vXIST_full_3



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -fm -pb -o v43_full_4

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_full_4.csv -o v43vXIST_full_4



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -fm -pb -o v43_full_5

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_full_5.csv -o v43vXIST_full_5



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -fm -pb -o v43_full_6

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_full_6.csv -o v43vXIST_full_6



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -fm -pb -o v43_full_7

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_full_7.csv -o v43vXIST_full_7



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -fm -pb -o v43_full_8

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_full_8.csv -o v43vXIST_full_8



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -fm -pb -o v43_full_9

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_full_9.csv -o v43vXIST_full_9



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -fm -pb -o v43_full_10

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_full_10.csv -o v43vXIST_full_10





