#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=2-0
#SBATCH --mem=96g



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 1000000 -fm -pb -o v43_1M_1

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_1M_1.csv -o v43vXIST_1M_1

seekr_adj_pval v43vXIST_1M_1.csv fdr_bh -o v43vXIST_1M_bh_1


seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 1000000 -fm -pb -o v43_1M_2

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_1M_2.csv -o v43vXIST_1M_2

seekr_adj_pval v43vXIST_1M_2.csv fdr_bh -o v43vXIST_1M_bh_2


seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 1000000 -fm -pb -o v43_1M_3

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_1M_3.csv -o v43vXIST_1M_3

seekr_adj_pval v43vXIST_1M_3.csv fdr_bh -o v43vXIST_1M_bh_3


seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 1000000 -fm -pb -o v43_1M_4

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_1M_4.csv -o v43vXIST_1M_4

seekr_adj_pval v43vXIST_1M_4.csv fdr_bh -o v43vXIST_1M_bh_4


seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 1000000 -fm -pb -o v43_1M_5

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_1M_5.csv -o v43vXIST_1M_5

seekr_adj_pval v43vXIST_1M_5.csv fdr_bh -o v43vXIST_1M_bh_5


seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 1000000 -fm -pb -o v43_1M_6

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_1M_6.csv -o v43vXIST_1M_6

seekr_adj_pval v43vXIST_1M_6.csv fdr_bh -o v43vXIST_1M_bh_6


seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 1000000 -fm -pb -o v43_1M_7

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_1M_7.csv -o v43vXIST_1M_7

seekr_adj_pval v43vXIST_1M_7.csv fdr_bh -o v43vXIST_1M_bh_7



seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 1000000 -fm -pb -o v43_1M_8

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_1M_8.csv -o v43vXIST_1M_8

seekr_adj_pval v43vXIST_1M_8.csv fdr_bh -o v43vXIST_1M_bh_8


seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 1000000 -fm -pb -o v43_1M_9

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_1M_9.csv -o v43vXIST_1M_9

seekr_adj_pval v43vXIST_1M_9.csv fdr_bh -o v43vXIST_1M_bh_9


seekr_find_dist v43_canonical.fa -k 6 -mdl common10 -sbt -sbs 1000000 -fm -pb -o v43_1M_10

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_1M_10.csv -o v43vXIST_1M_10

seekr_adj_pval v43vXIST_1M_10.csv fdr_bh -o v43vXIST_1M_bh_10




