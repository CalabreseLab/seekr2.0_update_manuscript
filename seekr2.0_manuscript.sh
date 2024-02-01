seekr_download_gencode lncRNA -g -r 43 -s human

seekr_filter_gencode v43_lncRNA.fa -gtf v43_lncRNA.chr_patch_hapl_scaff.annotation.gtf -len 500 -can -rd -o v43_canonical

seekr_norm_vectors v43_canonical.fa -k 6 -mv mean_6mers.npy -sv std_6mers.npy

seekr_kmer_counts XIST.fa -k 6 -o XIST_6mers.csv -mv mean_6mers.npy -sv std_6mers.npy

seekr_kmer_counts KCNQ1OT1.fa -k 6 -o OT1_6mers.csv -mv mean_6mers.npy -sv std_6mers.npy

seekr_pearson XIST_6mers.csv OT1_6mers.csv -o XIST_vs_OT1.csv

seekr_find_dist v43_canonical.fa -k 6 -sbt -sbs 100000 -fm -pb -pf modelfit6.pdf -o v43_can_6dists

seekr_find_dist v43_canonical.fa -k 6 -sbt -sbs 100000 -fm -pb -pf modelfit6.pdf -o v43_can_6dists

seekr_find_pval XIST.fa KCNQ1OT1.fa mean_6mers.npy std_6mers.npy 6 v43_can_6dists.csv -o XIST_v_OT1_pval

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_can_6dists.csv -o v43_v_XIST_pval

seekr_adj_pval v43_v_XIST_pval.csv fdr_bh -o v43_v_XIST_pvals_bh

seekr_find_dist v43_can500_namedchunks_500.fa -k 4 -sbt -sbs 100000 -fm -pb -pf modelfit4.pdf -o v43_can_4dists

seekr_norm_vectors v43_can500_namedchunks_500.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_find_pval XIST_manual_chunks.fa DLX6-AS1_namedchunks_500.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -ft distribution -bf 1 -o DLX6_v_X_pvals -pb

seekr_find_pval XIST_manual_chunks.fa LINC00632_namedchunks_500.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -ft distribution -bf 1 -o LINC00632_v_X_pvals -pb

seekr_find_pval XIST_manual_chunks.fa PCDH10-DT_namedchunks_500.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -ft distribution -bf 1 -o PCDH10-DT_v_X_pvals -pb

seekr_kmer_heatmap DLX6_v_X_pvals.csv 0 1 -th 0.05 -cr '#1b7837,#ffffff,#c51b7d' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o DLX6_v_X_pvals -hf pdf -hd 300

seekr_kmer_heatmap LINC00632_v_X_pvals.csv 0 1 -th 0.05 -cr '#1b7837,#ffffff,#c51b7d' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o LINC00632_v_X_pvals -hf pdf -hd 300

seekr_kmer_heatmap PCDH10-DT_v_X_pvals.csv 0 1 -th 0.05 -cr '#1b7837,#ffffff,#c51b7d' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o PCDH10-DT_v_X_pvals -hf pdf -hd 300

seekr_norm_vectors v43_canonical.fa -k 3 -mv mean_3mers.npy -sv std_3mers.npy

seekr_kmer_counts select_lncs.fa -k 3 -mv mean_3mers.npy -sv std_3mers.npy -o select_lnc_3mers.csv

seekr_kmer_heatmap select_lnc_3mers.csv 0 3 -th 1 -cr '#fcfc03,#ffffff,#1907e6' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o select_lnc_3mers -hf pdf -hd 300

seekr_norm_vectors v43_canonical.fa -k 3 -mv mean_3mers.npy -sv std_3mers.npy

seekr_kmer_msd_barplot select_lncs.fa mean_3mers.npy std_3mers.npy 3 -l Log2.post -ss mean -sm descending -tn 10 -xls 20 -yls 20 -xts 20 -yts 20 -o select_lncs_msd_barplot3 -pf pdf -d 300

seekr_norm_vectors v43_canonical.fa -k 3 -mv mean_3mers.npy -sv std_3mers.npy

seekr_kmer_count_barplot select_lncs.fa mean_3mers.npy std_3mers.npy 3 -l Log2.post -sm descending -tn 10 -xls 20 -yls 20 -xts 20 -yts 20 -ls 12 -o select_lncs_barplot3 -pf pdf -d 300

seekr_norm_vectors v43_can500_namedchunks_500.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_find_pval XIST_manual_chunks.fa XIST_manual_chunks.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -ft distribution -bf 1 -o X_v_X_pvals -pb

seekr_kmer_heatmap X_v_X_pvals.csv 0 1 -th 0.05 -cr '#1b7837,#ffffff,#c51b7d' -distm correlation -linkm complete -wratio 0.4 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o X_v_X_pvals -hf pdf -hd 300

seekr_kmer_dendrogram X_v_X_pvals.csv -dd column -distm correlation -linkm complete -ph 8 -wratio 0.5 -lfs 16 -o X_v_X_dendrogram -pf pdf -d 300

seekr_norm_vectors v43_can500_namedchunks_500.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_kmer_leiden XIST_manual_chunks.fa mean_4mers.npy std_4mers.npy 4 -a RBERVertexPartition -r 1.1 -pco 0.15 -cf XIST_manual_chunks 

seekr_norm_vectors v43_can500_namedchunks_500.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_find_pval v43_can500_namedchunks_500.fa XIST_repeats.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -ft distribution -bf 1 -o v43_chunks_v_ABDEF_pvals -pb


