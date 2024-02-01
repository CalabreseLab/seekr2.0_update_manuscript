# Seekr2.0 update manuscript

This github repo contains the codes needed to generate plots in the seekr 2.0 update manuscript.

## Installation

 To use this library, you have to have >=Python3.9.5 on your computer.

 Once you have Python, run:

 ```
 $ pip install seekr
 ```

 which will make both the command line tool and the python module available. Here listed are the codes to generate plots in console and in python separately.


## Console

### Figure 1

#### Panel A

For example, selecting the set of deduplicated “Ensembl_canonical” lncRNAs from GENCODE’s v43 lncRNA collection reduces the total number of lncRNA transcripts from 58,023 to 15,550.

```
seekr_download_gencode lncRNA -g -r 43 -s human

seekr_filter_gencode v43_lncRNA.fa -gtf v43_lncRNA.chr_patch_hapl_scaff.annotation.gtf -len 500 -can -rd -o v43_canonical
```

For example, comparing the human lncRNA XIST to that of another human lncRNA, KCNQ1OT1, using k-mer length of k=6 and the set of deduplicated “Ensembl_canonical” lncRNAs in GENCODE as a background set, we find that SEEKR returns a Pearson’s r value of 0.07, meaning that relative to the set of Ensembl_canonical lncRNAs the 6-mer profiles of XIST and KCNQ1OT1 exhibit a weak positive correlation (FIG).

```
seekr_norm_vectors v43_canonical.fa -k 6 -mv mean_6mers.npy -sv std_6mers.npy

seekr_kmer_counts XIST.fa -k 6 -o XIST_6mers.csv -mv mean_6mers.npy -sv std_6mers.npy

seekr_kmer_counts KCNQ1OT1.fa -k 6 -o OT1_6mers.csv -mv mean_6mers.npy -sv std_6mers.npy

seekr_pearson XIST_6mers.csv OT1_6mers.csv -o XIST_vs_OT1.csv
```

#### Panel B

SEEKR-derived Pearson’s r values for large sets of background sequences are often well-fit by log-normal distributions (Figure 1B).

```
seekr_find_dist v43_canonical.fa -k 6 -sbt -sbs 100000 -fm -pb -pf modelfit6.pdf -o v43_can_6dists
```

#### Panel C

Continuing from the example above, we find that XIST and KCNQ1OT1 are more similar to each other than ~93% of sampled pairwise comparisons of GENCODE canonical lncRNAs, despite their low positive Pearson’s r value (p-value of 0.07; Figure 1C).

Moreover, with two additional commands, we find some ~1500 GENCODE canonical lncRNAs whose overall k-mer contents are more similar to XIST than KCNQ1OT1 (Supp_table). 

The top ten most XIST-similar lncRNAs, along with their Benjamini-Hochberg adjusted p-values are shown in Table 1.

```
seekr_find_dist v43_canonical.fa -k 6 -sbt -sbs 100000 -fm -pb -pf modelfit6.pdf -o v43_can_6dists

seekr_find_pval XIST.fa KCNQ1OT1.fa mean_6mers.npy std_6mers.npy 6 v43_can_6dists.csv -o XIST_v_OT1_pval

seekr_find_pval v43_canonical.fa XIST.fa mean_6mers.npy std_6mers.npy 6 v43_can_6dists.csv -o v43_v_XIST_pval

seekr_adj_pval v43_v_XIST_pval.csv fdr_bh -o v43_v_XIST_pvals_bh
```

### Figure 2 graphing with chunks

#### generate chunks

Therefore, to determine whether our select XIST-like lncRNAs harbored notable regional similarities with XIST, we fragmented the lncRNAs into ~500 nucleotide (nt) chunks and compared the chunks in each XIST-like lncRNA to each chunk within XIST.

For XIST, we fragemented it based on its sequence features: each repeat is by itself an individual chunk. For the intervals that is longer than 1000nt, we fragmented it in the same way as the other lncRNAs.

The python code for this part is available in this github repository named: generate_chunks.py
The fragmented lncRNAs and XIST fasta file are also available: v43_can500_namedchunks_500.fa and XIST_manual_chunks.fa

#### Panel B to D

Because data suggest Repeats A, B, D, E, and F are discrete functional domains, we evaluated each XIST Repeat as a single intact chunk and separated the intervening XIST intervals into ~500 chunks. From these analyses, we found that LINC00632, DLX6-AS1, and PCDH10-DT each harbor significant similarity to Repeat E but lack similarity to the other Repeats in XIST (Figures 2B-D). Moreover, LINC00632, DLX6-AS1, and PCDH10-DT each harbor significant similarity to chunks distributed across the final exon of XIST (Figures 2B-D).

```
seekr_find_dist v43_can500_namedchunks_500.fa -k 4 -sbt -sbs 100000 -fm -pb -pf modelfit4.pdf -o v43_can_4dists

seekr_norm_vectors v43_can500_namedchunks_500.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_find_pval XIST_manual_chunks.fa DLX6-AS1_namedchunks_500.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -ft distribution -bf 1 -o DLX6_v_X_pvals -pb

seekr_find_pval XIST_manual_chunks.fa LINC00632_namedchunks_500.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -ft distribution -bf 1 -o LINC00632_v_X_pvals -pb

seekr_find_pval XIST_manual_chunks.fa PCDH10-DT_namedchunks_500.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -ft distribution -bf 1 -o PCDH10-DT_v_X_pvals -pb

seekr_kmer_heatmap DLX6_v_X_pvals.csv 0 1 -th 0.05 -cr '#1b7837,#ffffff,#c51b7d' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o DLX6_v_X_pvals -hf pdf -hd 300

seekr_kmer_heatmap LINC00632_v_X_pvals.csv 0 1 -th 0.05 -cr '#1b7837,#ffffff,#c51b7d' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o LINC00632_v_X_pvals -hf pdf -hd 300

seekr_kmer_heatmap PCDH10-DT_v_X_pvals.csv 0 1 -th 0.05 -cr '#1b7837,#ffffff,#c51b7d' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o PCDH10-DT_v_X_pvals -hf pdf -hd 300
```

#### Panel E to G

A closer examination revealed that these similarities could be attributed to the uniform enrichment of k-mers rich in A and T nucleotides; whereas, k-mers rich in G and C nucleotides were among the most variably enriched (Figures 2E-G).

Code for Panel E:

```
seekr_norm_vectors v43_canonical.fa -k 3 -mv mean_3mers.npy -sv std_3mers.npy

seekr_kmer_counts select_lncs.fa -k 3 -mv mean_3mers.npy -sv std_3mers.npy -o select_lnc_3mers.csv

seekr_kmer_heatmap select_lnc_3mers.csv 0 3 -th 1 -cr '#fcfc03,#ffffff,#1907e6' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o select_lnc_3mers -hf pdf -hd 300
```

Code for Panel F:

```
seekr_norm_vectors v43_canonical.fa -k 3 -mv mean_3mers.npy -sv std_3mers.npy

seekr_kmer_msd_barplot select_lncs.fa mean_3mers.npy std_3mers.npy 3 -l Log2.post -ss mean -sm descending -tn 10 -xls 20 -yls 20 -xts 20 -yts 20 -o select_lncs_msd_barplot3 -pf pdf -d 300
```

Code for Panel G:

```
seekr_norm_vectors v43_canonical.fa -k 3 -mv mean_3mers.npy -sv std_3mers.npy

seekr_kmer_count_barplot select_lncs.fa mean_3mers.npy std_3mers.npy 3 -l Log2.post -sm descending -tn 10 -xls 20 -yls 20 -xts 20 -yts 20 -ls 12 -o select_lncs_barplot3 -pf pdf -d 300
```

#### Panel H I and S1

We also unexpectedly observed that many of the 500nt chunks within XIST’s final exon were significantly more similar to each other than would be expected by chance (Figures 2H, 2I, and S1).

Code for Panel H:

```
seekr_norm_vectors v43_can500_namedchunks_500.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_find_pval XIST_manual_chunks.fa XIST_manual_chunks.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -ft distribution -bf 1 -o X_v_X_pvals -pb

seekr_kmer_heatmap X_v_X_pvals.csv 0 1 -th 0.05 -cr '#1b7837,#ffffff,#c51b7d' -distm correlation -linkm complete -wratio 0.4 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o X_v_X_pvals -hf pdf -hd 300
```

Code for Panel I:

```
seekr_kmer_dendrogram X_v_X_pvals.csv -dd column -distm correlation -linkm complete -ph 8 -wratio 0.5 -lfs 16 -o X_v_X_dendrogram -pf pdf -d 300
```

Code for Panel S1:

```
seekr_norm_vectors v43_can500_namedchunks_500.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_kmer_leiden XIST_manual_chunks.fa mean_4mers.npy std_4mers.npy 4 -a RBERVertexPartition -r 1.1 -pco 0.15 -cf XIST_manual_chunks 
```

The result files: XIST_manual_chunks_nodes_leiden.csv and XIST_manual_chunks_edges_leiden.csv are then directly imported into Gephi as nodes and edges files for network plotting with Force Atalas2 layout.

### Table 2 and Table S2

Given that the tandem repeats within XIST are some of the regions most essential for its ability to induce and maintain gene silencing, we used SEEKR to perform a parallel search for XIST-like lncRNAs among the set of GENCODE canonical lncRNAs. In this search, we separated all 15,993 GENCODE canonical lncRNAs into ~500nt fragments, and then used SEEKR to identify lncRNAs containing fragments that harbored significant k-mer similarity to XIST Repeats A, B, D, E, and F (p-value of <0.05).

We then summed the number of XIST-similar fragments in each lncRNA and used these sums to rank lncRNAs by their overall XIST-likeness. This analysis identified several intriguing lncRNAs (TABLE; sup table).

XIST_repeats.fa is also provided under this github repository.

```
seekr_norm_vectors v43_can500_namedchunks_500.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_find_pval v43_can500_namedchunks_500.fa XIST_repeats.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -ft distribution -bf 1 -o v43_chunks_v_ABDEF_pvals -pb
```

The result file v43_chunks_v_ABDEF_pvals.csv is then processed in R with the code: chunkABDE_organize.R to produce the count file: v43_lncRNA_ABDEF_count.csv. This file is then used to generate Table 2 and Table S2



















## Additional considerations

This section is a "not-so-quickstart", providing more complete views on selection of input data and parameter selection.

Some general advice for thinking about how to use SEEKR.
One challenge that we continually face in the lab is there are few ground truths in the lncRNA field and thus it is often unclear how to decide on the best parameters for sequence comparisons using SEEKR.
Below are some points that may be useful – these are also discussed in the conclusions of [PMID 31097619](https://pubmed.ncbi.nlm.nih.gov/31097619/)

### Selection of a set of sequences to use for the calculation of standardization vectors

In our experience, one of the most useful features of SEEKR is that it provides a metric of relative similarity.
Consider two lncRNAs, lncRNA-X, which is from mouse, and lncRNA-Y, which is from human.
One way to compare these two lncRNAs using SEEKR would be to calculate their *k*-mer profiles and compare these profiles in relation to all known mouse lncRNAs.
To do this, one would first use all mouse lncRNAs as an input for the "seekr_norm_vectors" script, to determine the mean and standard deviation of counts for all *k*-mers in all mouse lncRNAs.
Then, users would take those standardization vectors along with the sequences of lncRNA-X and lncRNA-Y and use them in the “kmer_counts” script to calculate *k*-mer profiles of lncRNA-X and lncRNA-Y in relation to all mouse lncRNAs.
Finally, users would employ the “seekr_pearson” script to determine how similar lncRNA-X and lncRNA-Y were to each other relative to all other mouse lncRNAs.
The point here is that the set of sequences used to calculate standardization vectors is a key variable – perhaps analogous to a reference gene in a quantitative PCR experiment or a loading control in a western blot.
Users should think through what comparison they are interested in performing and choose their set of standardization sequences accordingly.
In the example above, where all mouse lncRNAs were used for standardization, users might discover that “At the level of *k*-mers, lncRNA-X is more similar to lncRNA-Y than it is similar to 99% of other mouse lncRNAs”.
But changing the set of sequences for standardization can change the question being asked.
For example, users might be interested in comparing lncRNA-X and lncRNA-Y relative to all known human enhancer RNAs (eRNAs).
Here, using all human eRNAs as the standardization set, users might make an additional insightful discovery, perhaps: “At the level of *k*-mers, lncRNA-X is no more similar to lncRNA-Y than it is similar to the average human eRNA”.
Relatedly, when users are comparing two large groups of sequences (let us call these “set A” and “set B”), we again recommend thinking about what set of reference sequences would be best for standardization.
In most cases, users will probably want to use the same set of reference sequences to standardize *k*-mer counts in set A and in set B.  
Perhaps set A and set B should again be standardized relative to all mouse lncRNAs; or, both set A and set B be should be standardized relative to the sequences in set B.
But standardizing set A relative to itself, then standardizing set B relative to itself, then using “seekr_pearson” to compare the two sets of sequences might yield a non-sensical comparison.

### Selection of *k*-mer length

In our experience, the most robust biological trends have been relatively insensitive to the length of *k*-mer used in SEEKR.
Still, when deciding on a length of k to use for comparisons, we recommend using a *k*-mer length for which 4^k is similar to the length of the average feature or key feature that is being compared.
The reason for this is that as the length of k increases, so does the number of zero values for *k*-mer counts in a given sequence.
For example, there are 16384 possible 7-mers.
If users were interested in finding lncRNAs that are similar to lncRNA-X, which is 500 nucleotides long, a *k*-mer length of 7 would not be ideal, because the vector of 7-mer counts that corresponds to lncRNA-X would be dominated by zero values.
In this example, unless users had a specific rationale for searching 7-mers, a *k*-mer length of 4 (256 possible *k*-mers) or 5 (1024 possible *k*-mers) would provide the basis for a stronger comparison.

## Issues and Help

If you have questions about how you can use seekr in your own research, please send an email to jmcalabr@med.unc.edu

For full documentation, run:

```
$ seekr
```

Any suggestions, questions, or problems can be directed to our
[GitHub Issues page](https://github.com/CalabreseLab/seekr/issues).

Please also see [the pre-print](https://github.com/CalabreseLab/seekr/blob/logchanges/methods_mol_bio_seekr-v20.pdf) to a methods paper we wrote last year. This paper was originally scheduled to appear in Methods in Molecular Biology in 2020 but its publication date may be delayed.

## Citation

If you use this work, please cite:

```
Kirk, J. M., Kim, S. O., Inoue, K., Smola, M. J., Lee, D. M., Schertzer, M. D., … Calabrese, J. M. (2018). Functional classification of long non-coding RNAs by k -mer content. Nature Genetics, 50(10), 1474–1482. https://doi.org/10.1038/s41588-018-0207-8
```
