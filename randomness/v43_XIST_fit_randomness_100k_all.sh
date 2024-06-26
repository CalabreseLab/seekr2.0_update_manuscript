#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=2-0
#SBATCH --mem=96g



seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_1

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_1.csv -o v43vXIST_100k_1


seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_2

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_2.csv -o v43vXIST_100k_2


seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_3

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_3.csv -o v43vXIST_100k_3


seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_4

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_4.csv -o v43vXIST_100k_4


seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_5

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_5.csv -o v43vXIST_100k_5


seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_6

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_6.csv -o v43vXIST_100k_6


seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_7

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_7.csv -o v43vXIST_100k_7


seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_8

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_8.csv -o v43vXIST_100k_8


seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_9

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_9.csv -o v43vXIST_100k_9


seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_10

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_10.csv -o v43vXIST_100k_10



