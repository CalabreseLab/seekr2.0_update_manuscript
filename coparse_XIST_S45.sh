# XIST vs S4 and S5 

seekr_kmer_counts XIST.fa -k 6 -o XIST_6mers.csv -mv v25_mean6.npy -sv v25_std6.npy

seekr_kmer_counts v25_S4.fa -k 6 -o v25_S4_6mers.csv -mv v25_mean6.npy -sv v25_std6.npy

seekr_kmer_counts v25_S5.fa -k 6 -o v25_S5_6mers.csv -mv v25_mean6.npy -sv v25_std6.npy

seekr_kmer_counts vM10_S4.fa -k 6 -o vM10_S4_6mers.csv -mv v25_mean6.npy -sv v25_std6.npy

seekr_kmer_counts vM10_S5.fa -k 6 -o vM10_S5_6mers.csv -mv v25_mean6.npy -sv v25_std6.npy


seekr_pearson v25_S4_6mers.csv XIST_6mers.csv -o v25_S4_vs_XIST_rval.csv

seekr_pearson v25_S5_6mers.csv XIST_6mers.csv -o v25_S5_vs_XIST_rval.csv

seekr_pearson vM10_S4_6mers.csv XIST_6mers.csv -o vM10_S4_vs_XIST_rval.csv

seekr_pearson vM10_S5_6mers.csv XIST_6mers.csv -o vM10_S5_vs_XIST_rval.csv


seekr_find_pval v25_S4.fa XIST.fa v25_mean6.npy v25_std6.npy 6 v25_can_6dist.csv -o v25_S4_v_XIST_pval

seekr_adj_pval v25_S4_v_XIST_pval.csv fdr_bh -o v25_S4_v_XIST_pval_bh

seekr_find_pval v25_S5.fa XIST.fa v25_mean6.npy v25_std6.npy 6 v25_can_6dist.csv -o v25_S5_v_XIST_pval

seekr_adj_pval v25_S5_v_XIST_pval.csv fdr_bh -o v25_S5_v_XIST_pval_bh


seekr_find_pval vM10_S4.fa XIST.fa v25_mean6.npy v25_std6.npy 6 v25_can_6dist.csv -o vM10_S4_v_XIST_pval

seekr_adj_pval vM10_S4_v_XIST_pval.csv fdr_bh -o vM10_S4_v_XIST_pval_bh

seekr_find_pval vM10_S5.fa XIST.fa v25_mean6.npy v25_std6.npy 6 v25_can_6dist.csv -o vM10_S5_v_XIST_pval

seekr_adj_pval vM10_S5_v_XIST_pval.csv fdr_bh -o vM10_S5_v_XIST_pval_bh





# XIST repeats vs S4 and S5 chunks
seekr_kmer_counts XIST_repeats.fa -k 6 -o XIST_repeats_6mers.csv -mv v25_mean6.npy -sv v25_std6.npy

seekr_kmer_counts v25_S4_500chunk.fa -k 6 -o v25_S4_500chunk_6mers.csv -mv v25_mean6.npy -sv v25_std6.npy

seekr_kmer_counts v25_S5_500chunk.fa -k 6 -o v25_S5_500chunk_6mers.csv -mv v25_mean6.npy -sv v25_std6.npy

seekr_kmer_counts vM10_S4_500chunk.fa -k 6 -o vM10_S4_500chunk_6mers.csv -mv v25_mean6.npy -sv v25_std6.npy

seekr_kmer_counts vM10_S5_500chunk.fa -k 6 -o vM10_S5_500chunk_6mers.csv -mv v25_mean6.npy -sv v25_std6.npy


seekr_pearson v25_S4_500chunk_6mers.csv XIST_repeats_6mers.csv -o v25_S4_vs_XIST_repeats_rval.csv

seekr_pearson v25_S5_500chunk_6mers.csv XIST_repeats_6mers.csv -o v25_S5_vs_XIST_repeats_rval.csv

seekr_pearson vM10_S4_500chunk_6mers.csv XIST_repeats_6mers.csv -o vM10_S4_vs_XIST_repeats_rval.csv

seekr_pearson vM10_S5_500chunk_6mers.csv XIST_repeats_6mers.csv -o vM10_S5_vs_XIST_repeats_rval.csv


seekr_find_pval v25_S4_500chunk.fa XIST_repeats.fa v25_mean6.npy v25_std6.npy 6 v25_can_6dist.csv -o v25_S4_v_XIST_repeats_pval

seekr_adj_pval v25_S4_v_XIST_repeats_pval.csv fdr_bh -o v25_S4_v_XIST_repeats_pval_bh

seekr_find_pval v25_S5_500chunk.fa XIST_repeats.fa v25_mean6.npy v25_std6.npy 6 v25_can_6dist.csv -o v25_S5_v_XIST_repeats_pval

seekr_adj_pval v25_S5_v_XIST_repeats_pval.csv fdr_bh -o v25_S5_v_XIST_repeats_pval_bh


seekr_find_pval vM10_S4_500chunk.fa XIST_repeats.fa v25_mean6.npy v25_std6.npy 6 v25_can_6dist.csv -o vM10_S4_v_XIST_repeats_pval

seekr_adj_pval vM10_S4_v_XIST_repeats_pval.csv fdr_bh -o vM10_S4_v_XIST_repeats_pval_bh

seekr_find_pval vM10_S5_500chunk.fa XIST_repeats.fa v25_mean6.npy v25_std6.npy 6 v25_can_6dist.csv -o vM10_S5_v_XIST_repeats_pval

seekr_adj_pval vM10_S5_v_XIST_repeats_pval.csv fdr_bh -o vM10_S5_v_XIST_repeats_pval_bh




