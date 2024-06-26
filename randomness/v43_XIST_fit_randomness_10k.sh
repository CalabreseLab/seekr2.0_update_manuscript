#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=2-0
#SBATCH --mem=96g



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 10000 -fm -pb -o v43_10k_1

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_10k_1.csv -o v43vXIST_10k_1



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 10000 -fm -pb -o v43_10k_2

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_10k_2.csv -o v43vXIST_10k_2



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 10000 -fm -pb -o v43_10k_3

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_10k_3.csv -o v43vXIST_10k_3



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 10000 -fm -pb -o v43_10k_4

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_10k_4.csv -o v43vXIST_10k_4



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 10000 -fm -pb -o v43_10k_5

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_10k_5.csv -o v43vXIST_10k_5



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 10000 -fm -pb -o v43_10k_6

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_10k_6.csv -o v43vXIST_10k_6



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 10000 -fm -pb -o v43_10k_7

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_10k_7.csv -o v43vXIST_10k_7




seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 10000 -fm -pb -o v43_10k_8

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_10k_8.csv -o v43vXIST_10k_8



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 10000 -fm -pb -o v43_10k_9

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_10k_9.csv -o v43vXIST_10k_9



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 10000 -fm -pb -o v43_10k_10

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_10k_10.csv -o v43vXIST_10k_10





