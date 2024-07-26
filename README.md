# Seekr2.0 update manuscript

This github repo contains the codes needed to generate plots in the seekr 2.0 update manuscript.

## Installation

 To use this library, you have to have >=Python3.9.5 on your computer.

 Once you have Python, run:

 ```
 $ pip install seekr
 ```

 which will make both the command line tool and the python module available. Here listed are the codes to generate plots. Some plots are easier to generate with just the console commmands; while other are more complicated and involves using console commmands, python and R. 


## Figure 1

### Panel B and C

Select the set of deduplicated “Ensembl_canonical” lncRNAs that are greater than or equal to 500nt in length from the GENCODE v43 lncRNA collection.

```
seekr_download_gencode lncRNA -g -r 43 -s human

seekr_filter_gencode v43_lncRNA.fa -gtf v43_lncRNA.chr_patch_hapl_scaff.annotation.gtf -len 500 -can -rd -o v43_canonical
```

All codes for this part of analysis is saved under the **randomness folder**.

Use subset size = 100k, fit the data with 'all' models or 'common10' models. Calcualte the p values using the best fitted model. Codes saved as: **v43_XIST_fit_randomness_100k_all.sh** and **v43_XIST_fit_randomness_100k.sh**, respectively.

Compare results between 'all' best fitted model and 'common10' best fitted model. Generate the dataframe that includes the r value that corresponds to 95% percentile of each fitted model (where p=0.05), and the significant gene count using p value<0.05. Codes saved as **randomness_allvscommon10.py**. And plot the boxplot of the comparison results using **randomness_boxplot.R**.

Again, run v43_canonical.fa against XIST.fa at k=6 with different subset size (10k, 1M and all-no subsetting, 100k has been done in previous step) and fit the 'common10' models. Each simulate 10 times. Codes used here are saved as **v43_XIST_fit_randomness_10k/1M/full.sh**.

From the fitted model files, get the best fit model name and calculate the significant gene counts by p value <0.05 for each simulation. Plot the significante gene counts against different subset size using codes: **randomness_boxplot.R**.


### Panel D

Fit the data with k=4, 5 or 6 and check how kmer size affect the number of significant lncRNAs detected.

```
# k=4

seekr_norm_vectors v43_canonical.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_find_dist v43_canonical.fa -k 4 -sbt -sbs 100000 -fm -pb -pf modelfit4 -o v43_can_4dists

seekr_find_pval v43_canonical.fa XIST.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -o v43_v_XIST_pval_4

#k=5

seekr_norm_vectors v43_canonical.fa -k 5 -mv mean_5mers.npy -sv std_5mers.npy

seekr_find_dist v43_canonical.fa -k 5 -sbt -sbs 100000 -fm -pb -pf modelfit5 -o v43_can_5dists

seekr_find_pval v43_canonical.fa XIST.fa mean_5mers.npy std_5mers.npy 5 v43_can_5dists.csv -o v43_v_XIST_pval_5

# k=6

seekr_norm_vectors v43_canonical.fa -k 6 -mv mean_6mers.npy -sv std_6mers.npy

seekr_find_dist v43_canonical.fa -k 6 -sbt -sbs 100000 -fm -pb -pf modelfit6 -o v43_can_6dists

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_can_6dists.csv -o v43_v_XIST_pval_6

```
After getting the p values for k=4, 5 and 6, the boxplot was generated using **kmer_barplot.R**

### Panel E

The SEEKR console commands used to download GENCODE lncRNAs and calculate adjusted p-values for a sequence comparison of interest:

```
seekr_download_gencode lncRNA -g -r 43 -s human -fp v43_lncRNA.fa.gz -gp v43_lncRNA.gtf.gz

seekr_filter_gencode v43_lncRNA.fa -gtf v43_lncRNA.gtf -len 500 -can -rd -o v43_canonical

seekr_find_dist v43_canonical.fa -k 6 -sbt -fm -o v43_can_6dists

seekr_norm_vectors v43_canonical.fa -k 6 -mv mean_6mers.npy -sv std_6mers.npy

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_can_6dists.csv -o v43_v_XIST_pval

seekr_adj_pval v43_v_XIST_pval.csv fdr_bh -o v43_v_XIST_pvals_bh

```



### Table 1 and S1

Select the set of deduplicated “Ensembl_canonical” lncRNAs that are greater than or equal to 500nt in length from the GENCODE v43 lncRNA collection.

```
seekr_download_gencode lncRNA -g -r 43 -s human

seekr_filter_gencode v43_lncRNA.fa -gtf v43_lncRNA.chr_patch_hapl_scaff.annotation.gtf -len 500 -can -rd -o v43_canonical
```

Estimate a p value for the similarity between *XIST* and *KCNQ1OT1* using deduplicated GENCODE canonical lncRNAs as background and the best-fitting distribution from seekr_find_dist as the function to calculate significance.

Estimate p values for the similarity between *XIST* and all other GENCODE canonical lncRNAs.

```
seekr_norm_vectors v43_canonical.fa -k 6 -mv mean_6mers.npy -sv std_6mers.npy

seekr_find_dist v43_canonical.fa -k 6 -sbt -sbs 100000 -fm -pb -pf modelfit6 -o v43_can_6dists

seekr_find_pval XIST.fa KCNQ1OT1.fa mean_6mers.npy std_6mers.npy 6 v43_can_6dists.csv -o XIST_v_OT1_pval

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_can_6dists.csv -o v43_v_XIST_pval

seekr_adj_pval v43_v_XIST_pval.csv fdr_bh -o v43_v_XIST_pvals_bh
```

## Figure 2 graphing with chunks

### generate chunks

Therefore, to determine whether our select XIST-like lncRNAs harbored notable regional similarities with XIST, we fragmented the lncRNAs into ~500 nucleotide (nt) chunks and compared the chunks in each XIST-like lncRNA to each chunk within XIST.

For XIST, we fragemented it based on its sequence features: each repeat is by itself an individual chunk. For the intervals that is longer than 1000nt, we fragmented it in the same way as the other lncRNAs.

The python code for this part is available in this github repository named: **generate_chunks.py**

The fragmented lncRNAs and XIST fasta file are also available: **v43_can500_namedchunks_500.fa** and **XIST_manual_chunks.fa**

### Panel B to D

Fragment select *XIST*-like lncRNAs into ~500 nucleotide (nt) chunks and compare the chunks in each *XIST*-like lncRNA to each chunk within *XIST*. Evaluate each *XIST* Repeat as a single intact chunk and separate the intervening *XIST* intervals into ~500 chunks. 

```
seekr_find_dist v43_can500_namedchunks_500.fa -k 4 -sbt -sbs 100000 -fm -pb -pf modelfit4 -o v43_can_4dists

seekr_norm_vectors v43_can500_namedchunks_500.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_find_pval XIST_manual_chunks.fa DLX6-AS1_namedchunks_500.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -ft distribution -bf 1 -o DLX6_v_X_pvals -pb

seekr_find_pval XIST_manual_chunks.fa LINC00632_namedchunks_500.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -ft distribution -bf 1 -o LINC00632_v_X_pvals -pb

seekr_find_pval XIST_manual_chunks.fa PCDH10-DT_namedchunks_500.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -ft distribution -bf 1 -o PCDH10-DT_v_X_pvals -pb

seekr_kmer_heatmap DLX6_v_X_pvals.csv 0 1 -th 0.05 -cr '#1b7837,#ffffff,#c51b7d' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o DLX6_v_X_pvals -hf pdf -hd 300

seekr_kmer_heatmap LINC00632_v_X_pvals.csv 0 1 -th 0.05 -cr '#1b7837,#ffffff,#c51b7d' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o LINC00632_v_X_pvals -hf pdf -hd 300

seekr_kmer_heatmap PCDH10-DT_v_X_pvals.csv 0 1 -th 0.05 -cr '#1b7837,#ffffff,#c51b7d' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o PCDH10-DT_v_X_pvals -hf pdf -hd 300

```

### Panel E

Plot a heatmap of *k*-mer similarity scores across fragments of *XIST*-like lncRNAs of interest. 

```
seekr_download_gencode lncRNA -g -r 43 -s human

seekr_filter_gencode v43_lncRNA.fa -gtf v43_lncRNA.chr_patch_hapl_scaff.annotation.gtf -len 500 -can -rd -o v43_canonical

seekr_norm_vectors v43_canonical.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_kmer_counts Figure_2_select_3p_like_lncs_04_22_24.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy -o select_lnc_4mers.csv

seekr_kmer_heatmap select_lnc_4mers.csv 0 4 -th 1 -cr '#FCFC03,#FFFFFF,#1907E6' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o select_lnc_4mers -hf pdf -hd 300

```

### Panel F

Plot a barplot of most-enriched *k*-mers across fragments of *XIST*-like lncRNAs of interest. 

```
seekr_download_gencode lncRNA -g -r 43 -s human

seekr_filter_gencode v43_lncRNA.fa -gtf v43_lncRNA.chr_patch_hapl_scaff.annotation.gtf -len 500 -can -rd -o v43_canonical

seekr_norm_vectors v43_canonical.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_kmer_msd_barplot Figure_2_select_3p_like_lncs_04_22_24.fa mean_4mers.npy std_4mers.npy 4 -l Log2.post -ss mean -sm descending -tn 10 -xls 20 -yls 20 -xts 20 -yts 20 -o select_lncs_msd_barplot4 -pf pdf -d 300
```

### Panel G

Plot a barplot of most variable *k*-mers across fragments of *XIST*-like lncRNAs of interest. 

```
seekr_download_gencode lncRNA -g -r 43 -s human

seekr_filter_gencode v43_lncRNA.fa -gtf v43_lncRNA.chr_patch_hapl_scaff.annotation.gtf -len 500 -can -rd -o v43_canonical

seekr_norm_vectors v43_canonical.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_kmer_count_barplot Figure_2_select_3p_like_lncs_04_22_24.fa mean_4mers.npy std_4mers.npy 4 -l Log2.post -sm descending -tn 10 -xls 20 -yls 20 -xts 20 -yts 20 -ls 12 -o select_lncs_barplot4 -pf pdf -d 300

```

### Panel H

Compare *k*-mer similarity among fragments of *XIST* at *k*-mer length *k* = 4. Calculate p-values and plot them in a heatmap.

```
seekr_norm_vectors v43_can500_namedchunks_500.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_find_pval XIST_manual_chunks.fa XIST_manual_chunks.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -ft distribution -bf 1 -o X_v_X_pvals -pb

seekr_kmer_heatmap X_v_X_pvals.csv 0 1 -th 0.05 -cr '#1b7837,#ffffff,#c51b7d' -distm correlation -linkm complete -wratio 0.4 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o X_v_X_pvals -hf pdf -hd 300
```

### Panel I

Compare *k*-mer similarity among fragments of *XIST* at *k*-mer length *k* = 4. From p-values calculated above, create a dendrogram of fragment similarities across *XIST*.

```
seekr_kmer_dendrogram X_v_X_pvals.csv -dd column -distm correlation -linkm complete -ph 8 -wratio 0.5 -lfs 16 -o X_v_X_dendrogram -pf pdf -d 300
```

## Figure 4 human (s)eCLIP analysis

All codes for this section are saved under the **eclips folder**.

Get the values from Table_S2_v43_lncRNA_ABDEF_count_expression_2_19_2024_formatted.xlsx and calculate the mean (k_tot_mean) of k_tot_1 and k_tot_2, and filter for the genes that have k_tot_mean > 0.0625. Use this as the genelist to filter for chunks list in v43_chunks_v_ABDEF_pvals.csv (previous results) for chunks that belongs to a gene within the genelist. For each repeat in XIST (rA, rF, rB1, rB2, rD and rE), list the chunks that are significant for that specific repeat. 

For each of the significant chunk, use the chunkID to match it with a transcript in the corresponding gtf file (gencode.v43.long_noncoding_RNAs.gtf). First, locate the chunk along the transcript based on the coordinates in the chunkID (for example from 500 to 1000 of a transcript). Then get the exons that this transcript coordinates (500-1000) overlap. Lastly convert the transcript coordinates of the chunk to genomic coordinates, based on the relative transcript coordinates (500-1000) of the chunk, genomic coordinates of the overlapped exons and exon lengths. Then save it as a bedfile. If the chunk overlaps with multiple exons, it will be converted to multiple entries in the bedfile. In this way, a bedfile contains genomic coordinates of all significant chunks for each XIST repeat will be generated.

Randomization of these bedfiles is performed among the region covered by the genelist that has k_tot_mean > 0.0625. Each gene is then converted to a dataframe contains the transcriptID, geneID, name and its length (the pool). For each XIST repeat, the unique transcript IDs of all significant chunks are listed as the sigIDList. For each significant chunk, randomly choose an entry in the pool that has length greater than the chunk length and its transcript ID does not exist in the sigIDList (do not resample the genes that have a significant chunk). Then randomly choose the start coordinate from 0 to (transcript length-chunk length). Therefore for each randomly choosed region, its transcript ID, start coordinate (on the transcript level) and length is set.

After this transcript level randomization is done for each repeat of XIST, it is then converted to genomic coordinates using the same method mentioned above. In this way, a randomized control bedfile is generated for all significant chunks of each XIST repeat. 

All these are performed by the code: **eclips_bedfile_rand.R**.

Human (s)eCLIP bam files for RBM15, HNRNPK, MATR3, PTBP1 and HNRNPM are downloaded from [ENCODE](https://www.encodeproject.org/). HNRNPK are seCLIP (single end) while RBM15, PTBP1, HNRNPM, MATR3 are eCLIP (paired end). For bedtools multicov counting, bam file for HNRNPK should be aligned to the same strand (using argument -s). For RBM15, PTBP1, HNRNPM, MATR3, we can take only the second in pair and align to the same strand (using argument -s). In this way all bam files can use the same argument for multicov and can run altogether in the same command. Downloaded bam files are firstly sorted, indexed. Replicates are then merged for each RBP. For RBM15, PTBP1, HNRNPM, MATR3, extract the second read in pair and sort then index. In this way the bam files are prepared for bedtools multicov counting. 

Codes for this part is saved as: **eclips_sort_index_merge_bam.sh**.


Bedtools multicov counting is performed for all bedfiles of significant chunks for each XIST repeat and their corresponding randomized control bedfiles, using code: **eclips_multicov_counts_eCLIPs_chunk.sh**.


In the previous step if a chunk overlaps with multiple exons, it would be split into multiple entries in the bedfile according to the number of exons overlapped. Therefore for the output of multicov, counts for these split entries needs to be merged back to that chunk
so that the significant chunk counts could be compared with its randomized control. Then for each repeat and each RBP, paired Wilcoxon test (or Wilcoxon signed-rank test) is performed to determine whether the counts in the significant chunk is siginificant more than the counts in their corresponding control. The number of counts for each significant chunk and randomized control for each RBP and each repeat are plotted as a boxplot. Codes fo this part is saved as: **eclips_multicov_analysis_chunk.R**.


## Figure 6 coPARSE analysis

All codes for this part of analysis is saved under the **coparse folder**.

From [GENCODE](https://www.gencodegenes.org/) download lncRNA fasta and gtf file for human: v25 and mouse: vM10.

Filter both the human and mouse sequences and setting the background distribution of seekr r values using filtered human lncRNAs. Exract human and mouse sequences corresponds to Table S4 and Table S5a in the [coPARSE paper](https://www.nature.com/articles/s41588-023-01620-7). Calculate the seekr r values and p values for all human vs mouse pair in Table S4 and S5a. Codes for this part are saved as **coparse_analysis.py**. Here we did not use the command line function as there are functions used that are outside the seekr package and are not command line compatible. So it is easier to accomplish the whole session in Python. 

Plot seekr r value background distribution as a density plot; overlap with the human vs mouse r values in Table S4 and Table S5a as histogram. The histogram is then colored according to their p values. Codes for this part are saved as **coparse_plotting.R**. 

Cut human and mouse sequences in Table S4 and Table S5a into chunks, using code: **coparse_generate_chunks.py**.

Calculate the seekr r values and p values for sequences in Table S4 and S5a vs XIST full sequence; and chunked sequences of Table S4 and Table S5 vs XIST repeats. All functions needed here are within the seekr package. So here we provide the codes as in command line: **coparse_XIST_S45.sh**.

Analyze the results to get the significance chunk count per XIST repeat per gene and then combine all results together. Codes are saved as **coparse_chunk_analysis.R**.






## Figure S1

Regional similarities between XIST and LINC00632, DLX6-AS1, and PCDH10-DT at k-mer lengths k = 4, 5, and 6. Please refer to codes for Figure 2 for k =4.

```
seekr_download_gencode lncRNA -g -r 43 -s human

seekr_filter_gencode v43_lncRNA.fa -gtf v43_lncRNA.chr_patch_hapl_scaff.annotation.gtf -len 500 -can -rd -o v43_canonical

# k=5

seekr_norm_vectors v43_canonical.fa -k 5 -mv mean_5mers.npy -sv std_5mers.npy

seekr_find_pval ../XIST_manual_chunks_withF.fa ../DLX6-AS1_namedchunks_500.fa mean_5mers.npy std_5mers.npy 5 v43_can_5dists.csv -ft distribution -bf 1 -o DLX6_v_X_pvals -pb

seekr_find_pval ../XIST_manual_chunks_withF.fa ../LINC00632_namedchunks_500.fa mean_5mers.npy std_5mers.npy 5 v43_can_5dists.csv -ft distribution -bf 1 -o LINC00632_v_X_pvals -pb

seekr_find_pval ../XIST_manual_chunks_withF.fa ../PCDH10-DT_namedchunks_500.fa mean_5mers.npy std_5mers.npy 5 v43_can_5dists.csv -ft distribution -bf 1 -o PCDH10-DT_v_X_pvals -pb

seekr_kmer_heatmap DLX6_v_X_pvals.csv 0 1 -th 0.05 -cr '#1B7837,#FFFFFF,#C51B7D' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o DLX6_v_X_pvals -hf pdf -hd 300

seekr_kmer_heatmap LINC00632_v_X_pvals.csv 0 1 -th 0.05 -cr '#1B7837,#FFFFFF,#C51B7D' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o LINC00632_v_X_pvals -hf pdf -hd 300

seekr_kmer_heatmap PCDH10-DT_v_X_pvals.csv 0 1 -th 0.05 -cr '#1B7837,#FFFFFF,#C51B7D' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o PCDH10-DT_v_X_pvals -hf pdf -hd 300

# k=6

seekr_norm_vectors v43_canonical.fa -k 6 -mv mean_6mers.npy -sv std_6mers.npy

seekr_find_pval ../XIST_manual_chunks_withF.fa ../DLX6-AS1_namedchunks_500.fa mean_6mers.npy std_6mers.npy 6 v43_can_6dists.csv -ft distribution -bf 1 -o DLX6_v_X_pvals -pb

seekr_find_pval ../XIST_manual_chunks_withF.fa ../LINC00632_namedchunks_500.fa mean_6mers.npy std_6mers.npy 6 v43_can_6dists.csv -ft distribution -bf 1 -o LINC00632_v_X_pvals -pb

seekr_find_pval ../XIST_manual_chunks_withF.fa ../PCDH10-DT_namedchunks_500.fa mean_6mers.npy std_6mers.npy 6 v43_can_6dists.csv -ft distribution -bf 1 -o PCDH10-DT_v_X_pvals -pb

seekr_kmer_heatmap DLX6_v_X_pvals.csv 0 1 -th 0.05 -cr '#1B7837,#FFFFFF,#C51B7D' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o DLX6_v_X_pvals -hf pdf -hd 300

seekr_kmer_heatmap LINC00632_v_X_pvals.csv 0 1 -th 0.05 -cr '#1B7837,#FFFFFF,#C51B7D' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o LINC00632_v_X_pvals -hf pdf -hd 300

seekr_kmer_heatmap PCDH10-DT_v_X_pvals.csv 0 1 -th 0.05 -cr '#1B7837,#FFFFFF,#C51B7D' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o PCDH10-DT_v_X_pvals -hf pdf -hd 300
```

## Fig S2

Compare *k*-mer similarity among fragments of *XIST* at *k*-mer length *k* = 4 and use Leiden to assign fragments to communities.

```
seekr_norm_vectors v43_can500_namedchunks_500.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_kmer_leiden XIST_manual_chunks.fa mean_4mers.npy std_4mers.npy 4 -a RBERVertexPartition -r 1.1 -pco 0.15 -cf XIST_manual_chunks 
```

The result files: **XIST_manual_chunks_nodes_leiden.csv** and **XIST_manual_chunks_edges_leiden.csv** are then directly imported into Gephi as nodes and edges files for network plotting with Force Atalas2 layout.

## Table 2 and Table S2

Fragment GENCODE canonical lncRNAs into ~500nt chunks and compare each GENCODE lncRNA chunk to each *XIST* Repeat, using *k*-mer *k* = 4. Calculate p values of similarity using GENCODE lncRNA chunks as background and fitting to a lognormal distribution.

Define similarity as any fragment whose SEEKR-derived similarity to an *XIST* Repeat has a p value of less than 0.05. Tally the number of *XIST*-similar fragments in each lncRNA and used these sums to rank lncRNAs by their overall *XIST*-likeness.

**XIST_repeats.fa** is also provided under this github repository.

```
seekr_norm_vectors v43_can500_namedchunks_500.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_find_pval v43_can500_namedchunks_500.fa XIST_repeats.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -ft distribution -bf 1 -o v43_chunks_v_ABDEF_pvals -pb
```

The result file **v43_chunks_v_ABDEF_pvals.csv** is then processed in R with the code: **chunkABDE_organize.R** to produce the count file: **v43_lncRNA_ABDEF_count.csv**. This file is then used to generate Table 2 and Table S2





## Issues and Help

If you have questions about how you can use seekr in your own research, please send an email to jmcalabr@med.unc.edu

For full documentation, run:

```
$ seekr
```

For further details about the seekr package please refer to the
[Seekr GitHub page](https://github.com/CalabreseLab/seekr).

Please also see [the pre-print](https://github.com/CalabreseLab/seekr/blob/logchanges/methods_mol_bio_seekr-v20.pdf) to a methods paper we wrote last year. This paper was originally scheduled to appear in Methods in Molecular Biology in 2020 but its publication date may be delayed.

## Citation

If you use this work, please cite:

```
Kirk, J. M., Kim, S. O., Inoue, K., Smola, M. J., Lee, D. M., Schertzer, M. D., … Calabrese, J. M. (2018). Functional classification of long non-coding RNAs by k -mer content. Nature Genetics, 50(10), 1474–1482. https://doi.org/10.1038/s41588-018-0207-8
```

