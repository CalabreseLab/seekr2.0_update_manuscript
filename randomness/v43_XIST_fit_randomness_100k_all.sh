#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=2-0
#SBATCH --mem=96g



seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_1

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_1.csv -o v43vXIST_100k_1

seekr_adj_pval v43vXIST_100k_1.csv fdr_bh -o v43vXIST_100k_bh_1

# fit using the 2nd and 3rd bestfit
seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_1.csv -bf 2 -o v43vXIST_2_100k_1

seekr_adj_pval v43vXIST_2_100k_1.csv fdr_bh -o v43vXIST_2_100k_bh_1


seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_1.csv -bf 3 -o v43vXIST_3_100k_1

seekr_adj_pval v43vXIST_3_100k_1.csv fdr_bh -o v43vXIST_3_100k_bh_1

################

seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_2

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_2.csv -o v43vXIST_100k_2

seekr_adj_pval v43vXIST_100k_2.csv fdr_bh -o v43vXIST_100k_bh_2


seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_2.csv -bf 2 -o v43vXIST_2_100k_2

seekr_adj_pval v43vXIST_2_100k_2.csv fdr_bh -o v43vXIST_2_100k_bh_2

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_2.csv -bf 3 -o v43vXIST_3_100k_2

seekr_adj_pval v43vXIST_3_100k_2.csv fdr_bh -o v43vXIST_3_100k_bh_2

###################


seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_3

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_3.csv -o v43vXIST_100k_3

seekr_adj_pval v43vXIST_100k_3.csv fdr_bh -o v43vXIST_100k_bh_3


seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_3.csv -bf 2 -o v43vXIST_2_100k_3

seekr_adj_pval v43vXIST_2_100k_3.csv fdr_bh -o v43vXIST_2_100k_bh_3


seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_3.csv -bf 3 -o v43vXIST_3_100k_3

seekr_adj_pval v43vXIST_3_100k_3.csv fdr_bh -o v43vXIST_3_100k_bh_3



#############################
seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_4

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_4.csv -o v43vXIST_100k_4

seekr_adj_pval v43vXIST_100k_4.csv fdr_bh -o v43vXIST_100k_bh_4


seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_4.csv -bf 2 -o v43vXIST_2_100k_4

seekr_adj_pval v43vXIST_2_100k_4.csv fdr_bh -o v43vXIST_2_100k_bh_4


seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_4.csv -bf 3 -o v43vXIST_3_100k_4

seekr_adj_pval v43vXIST_3_100k_4.csv fdr_bh -o v43vXIST_3_100k_bh_4


####################################

seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_5

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_5.csv -o v43vXIST_100k_5

seekr_adj_pval v43vXIST_100k_5.csv fdr_bh -o v43vXIST_100k_bh_5


seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_5.csv -bf 2 -o v43vXIST_2_100k_5

seekr_adj_pval v43vXIST_2_100k_5.csv fdr_bh -o v43vXIST_2_100k_bh_5

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_5.csv -bf 3 -o v43vXIST_3_100k_5

seekr_adj_pval v43vXIST_3_100k_5.csv fdr_bh -o v43vXIST_3_100k_bh_5

########################
seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_6

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_6.csv -o v43vXIST_100k_6

seekr_adj_pval v43vXIST_100k_6.csv fdr_bh -o v43vXIST_100k_bh_6


seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_6.csv -bf 2 -o v43vXIST_2_100k_6

seekr_adj_pval v43vXIST_2_100k_6.csv fdr_bh -o v43vXIST_2_100k_bh_6

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_6.csv -bf 3 -o v43vXIST_3_100k_6

seekr_adj_pval v43vXIST_3_100k_6.csv fdr_bh -o v43vXIST_3_100k_bh_6


########################
seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_7

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_7.csv -o v43vXIST_100k_7

seekr_adj_pval v43vXIST_100k_7.csv fdr_bh -o v43vXIST_100k_bh_7


seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_7.csv -bf 2 -o v43vXIST_2_100k_7

seekr_adj_pval v43vXIST_2_100k_7.csv fdr_bh -o v43vXIST_2_100k_bh_7

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_7.csv -bf 3 -o v43vXIST_3_100k_7

seekr_adj_pval v43vXIST_3_100k_7.csv fdr_bh -o v43vXIST_3_100k_bh_7


########################
seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_8

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_8.csv -o v43vXIST_100k_8

seekr_adj_pval v43vXIST_100k_8.csv fdr_bh -o v43vXIST_100k_bh_8


seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_8.csv -bf 2 -o v43vXIST_2_100k_8

seekr_adj_pval v43vXIST_2_100k_8.csv fdr_bh -o v43vXIST_2_100k_bh_8


seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_8.csv -bf 3 -o v43vXIST_3_100k_8

seekr_adj_pval v43vXIST_3_100k_8.csv fdr_bh -o v43vXIST_3_100k_bh_8

########################
seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_9

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_9.csv -o v43vXIST_100k_9

seekr_adj_pval v43vXIST_100k_9.csv fdr_bh -o v43vXIST_100k_bh_9


seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_9.csv -bf 2 -o v43vXIST_2_100k_9

seekr_adj_pval v43vXIST_2_100k_9.csv fdr_bh -o v43vXIST_2_100k_bh_9


seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_9.csv -bf 3 -o v43vXIST_3_100k_9

seekr_adj_pval v43vXIST_3_100k_9.csv fdr_bh -o v43vXIST_3_100k_bh_9

########################
seekr_find_dist v43_canonical.fa -k 6 -mdl all -sbt -sbs 100000 -fm -pb -o v43_100k_10

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_10.csv -o v43vXIST_100k_10

seekr_adj_pval v43vXIST_100k_10.csv fdr_bh -o v43vXIST_100k_bh_10


seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_10.csv -bf 2 -o v43vXIST_2_100k_10

seekr_adj_pval v43vXIST_2_100k_10.csv fdr_bh -o v43vXIST_2_100k_bh_10

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_100k_10.csv -bf 3 -o v43vXIST_3_100k_10

seekr_adj_pval v43vXIST_3_100k_10.csv fdr_bh -o v43vXIST_3_100k_bh_10




