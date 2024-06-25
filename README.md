# Seekr2.0 update manuscript

This github repo contains the codes needed to generate plots in the seekr 2.0 update manuscript.

## Installation

 To use this library, you have to have >=Python3.9.5 on your computer.

 Once you have Python, run:

 ```
 $ pip install seekr
 ```

 which will make both the command line tool and the python module available. Here listed are the codes to generate plots in console and in python separately. seekr2.0_manuscript.sh and seekr2.0_manuscript.py contains the codes all togher without the comments.


## Console

### Figure 1

#### Panel A

Select the set of deduplicated “Ensembl_canonical” lncRNAs that are greater than or equal to 500nt in length from the GENCODE v43 lncRNA collection.

```
seekr_download_gencode lncRNA -g -r 43 -s human

seekr_filter_gencode v43_lncRNA.fa -gtf v43_lncRNA.chr_patch_hapl_scaff.annotation.gtf -len 500 -can -rd -o v43_canonical
```

Compare the human lncRNA *XIST* to that of another human lncRNA, *KCNQ1OT1*, using *k*-mer length of *k* = 6 and the set of deduplicated $\geq$ 500nt “Ensembl_canonical” lncRNAs in GENCODE as a background set.

```
seekr_norm_vectors v43_canonical.fa -k 6 -mv mean_6mers.npy -sv std_6mers.npy

seekr_kmer_counts XIST.fa -k 6 -o XIST_6mers.csv -mv mean_6mers.npy -sv std_6mers.npy

seekr_kmer_counts KCNQ1OT1.fa -k 6 -o OT1_6mers.csv -mv mean_6mers.npy -sv std_6mers.npy

seekr_pearson XIST_6mers.csv OT1_6mers.csv -o XIST_vs_OT1.csv

seekr_kmer_counts v43_canonical.fa -k 6 -o v43_6mers.csv -mv mean_6mers.npy -sv std_6mers.npy

seekr_pearson v43_6mers.csv XIST_6mers.csv -o v43_vs_XIST_rval.csv
```

#### Panel B

Fit Python’s ‘common10’ distributions to 100000 randomly selected pairwise comparisons of sequences from v43_canonical.fa using *k*-mer length *k* = 6.

```
seekr_find_dist v43_canonical.fa -k 6 -sbt -sbs 100000 -fm -pb -pf modelfit6.pdf -o v43_can_6dists
```

#### Panel C, Table 1 and S1

Estimate a p value for the similarity between *XIST* and *KCNQ1OT1* using deduplicated GENCODE canonical lncRNAs as background and the best-fitting distribution from seekr_find_dist as the function to calculate significance.


Estimate p values for the similarity between *XIST* and all other GENCODE canonical lncRNAs.


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

Fragment select *XIST*-like lncRNAs into ~500 nucleotide (nt) chunks and compare the chunks in each *XIST*-like lncRNA to each chunk within *XIST*. Evaluate each *XIST* Repeat as a single intact chunk and separate the intervening *XIST* intervals into ~500 chunks. 

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

#### Panel E

Plot a heatmap of *k*-mer similarity scores across fragments of *XIST*-like lncRNAs of interest. Use *k* = 3 for ease of visualizing individual *k*-mers.

```
seekr_norm_vectors v43_canonical.fa -k 3 -mv mean_3mers.npy -sv std_3mers.npy

seekr_kmer_counts select_lncs.fa -k 3 -mv mean_3mers.npy -sv std_3mers.npy -o select_lnc_3mers.csv

seekr_kmer_heatmap select_lnc_3mers.csv 0 3 -th 1 -cr '#fcfc03,#ffffff,#1907e6' -distm correlation -linkm complete -wratio 0.3 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o select_lnc_3mers -hf pdf -hd 300
```

#### Panel F

Plot a barplot of most-enriched *k*-mers across fragments of *XIST*-like lncRNAs of interest. Use *k* = 3 for ease of visualizing individual *k*-mers.

```
seekr_norm_vectors v43_canonical.fa -k 3 -mv mean_3mers.npy -sv std_3mers.npy

seekr_kmer_msd_barplot select_lncs.fa mean_3mers.npy std_3mers.npy 3 -l Log2.post -ss mean -sm descending -tn 10 -xls 20 -yls 20 -xts 20 -yts 20 -o select_lncs_msd_barplot3 -pf pdf -d 300
```

#### Panel G

Plot a barplot of most variable *k*-mers across fragments of *XIST*-like lncRNAs of interest. Use *k* = 3 for ease of visualizing individual *k*-mers.

```
seekr_norm_vectors v43_canonical.fa -k 3 -mv mean_3mers.npy -sv std_3mers.npy

seekr_kmer_count_barplot select_lncs.fa mean_3mers.npy std_3mers.npy 3 -l Log2.post -sm descending -tn 10 -xls 20 -yls 20 -xts 20 -yts 20 -ls 12 -o select_lncs_barplot3 -pf pdf -d 300
```

#### Panel H

Compare *k*-mer similarity among fragments of *XIST* at *k*-mer length *k* = 4. Calculate p-values and plot them in a heatmap.

```
seekr_norm_vectors v43_can500_namedchunks_500.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_find_pval XIST_manual_chunks.fa XIST_manual_chunks.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -ft distribution -bf 1 -o X_v_X_pvals -pb

seekr_kmer_heatmap X_v_X_pvals.csv 0 1 -th 0.05 -cr '#1b7837,#ffffff,#c51b7d' -distm correlation -linkm complete -wratio 0.4 -hratio 0.3 -xts 16 -yts 16 -cfs 16 -o X_v_X_pvals -hf pdf -hd 300
```

#### Panel I

Compare *k*-mer similarity among fragments of *XIST* at *k*-mer length *k* = 4. From p-values calculated above, create a dendrogram of fragment similarities across *XIST*.

```
seekr_kmer_dendrogram X_v_X_pvals.csv -dd column -distm correlation -linkm complete -ph 8 -wratio 0.5 -lfs 16 -o X_v_X_dendrogram -pf pdf -d 300
```

#### Panel S1

Compare *k*-mer similarity among fragments of *XIST* at *k*-mer length *k* = 4 and use Leiden to assign fragments to communities.

```
seekr_norm_vectors v43_can500_namedchunks_500.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_kmer_leiden XIST_manual_chunks.fa mean_4mers.npy std_4mers.npy 4 -a RBERVertexPartition -r 1.1 -pco 0.15 -cf XIST_manual_chunks 
```

The result files: XIST_manual_chunks_nodes_leiden.csv and XIST_manual_chunks_edges_leiden.csv are then directly imported into Gephi as nodes and edges files for network plotting with Force Atalas2 layout.

### Table 2 and Table S2

Fragment GENCODE canonical lncRNAs into ~500nt chunks and compare each GENCODE lncRNA chunk to each *XIST* Repeat, using *k*-mer *k* = 4. Calculate p values of similarity using GENCODE lncRNA chunks as background and fitting to a lognormal distribution.

Define similarity as any fragment whose SEEKR-derived similarity to an *XIST* Repeat has a p value of less than 0.05. Tally the number of *XIST*-similar fragments in each lncRNA and used these sums to rank lncRNAs by their overall *XIST*-likeness.

XIST_repeats.fa is also provided under this github repository.

```
seekr_norm_vectors v43_can500_namedchunks_500.fa -k 4 -mv mean_4mers.npy -sv std_4mers.npy

seekr_find_pval v43_can500_namedchunks_500.fa XIST_repeats.fa mean_4mers.npy std_4mers.npy 4 v43_can_4dists.csv -ft distribution -bf 1 -o v43_chunks_v_ABDEF_pvals -pb
```

The result file v43_chunks_v_ABDEF_pvals.csv is then processed in R with the code: chunkABDE_organize.R to produce the count file: v43_lncRNA_ABDEF_count.csv. This file is then used to generate Table 2 and Table S2



## Python

### Figure 1

#### Panel A

Select the set of deduplicated “Ensembl_canonical” lncRNAs that are greater than or equal to 500nt in length from the GENCODE v43 lncRNA collection.

```python
from seekr import fasta
from seekr import filter_gencode

downloader = fasta.Downloader()
downloader.get_gencode(biotype='lncRNA', species='human', gtf=True, release='43', unzip=True)

headers, seqs = filter_gencode.filter_gencode(fasta_path='v43_lncRNA.fa', 
                                              gtf_path='v43_lncRNA.chr_patch_hapl_scaff.annotation.gtf',
                                              len_threshold=500, canonical=True,
                                              rm_dup=True, outputname='v43_canonical')
```

Compare the human lncRNA *XIST* to that of another human lncRNA, *KCNQ1OT1*, using *k*-mer length of *k* = 6 and the set of deduplicated >=500nt “Ensembl_canonical” lncRNAs in GENCODE as a background set.

```python
from seekr.kmer_counts import BasicCounter
from seekr.pearson import pearson as seekrPearson
import numpy as np

bkg_norm_6 = BasicCounter('v43_canonical.fa', k=6)
bkg_norm_6.get_counts()

mean_path = 'mean_6mers.npy'
std_path = 'std_6mers.npy'
np.save(mean_path, bkg_norm_6.mean)
np.save(std_path, bkg_norm_6.std)


xist_count = BasicCounter(infasta='XIST.fa', outfile='XIST_6mers.csv', k=6,
                          mean='mean_6mers.npy', std='std_6mers.npy', 
                          log2='Log2.post')
xist_count. make_count_file()


ot1_count = BasicCounter(infasta='KCNQ1OT1.fa', outfile='OT1_6mers.csv', k=6, 
                         mean='mean_6mers.npy', std='std_6mers.npy', 
                         log2='Log2.post')
ot1_count. make_count_file()


sim = seekrPearson(xist_count.counts,ot1_count.counts, outfile='XIST_vs_OT1.csv')
```

#### Panel B

Fit Python’s ‘common10’ distributions to 100000 randomly selected pairwise comparisons of sequences from v43_canonical.fa using *k*-mer length *k* = 6.

```python
from seekr import find_dist

v43_6dists = find_dist.find_dist(inputseq='v43_canonical.fa', k_mer=6, 
                                 log2='Log2.post', models='common10', 
                                 subsetting=True, subset_size = 100000, 
                                 fit_model=True, statsmethod='ks',
                                 progress_bar=True, plotfit='modelfit6.pdf', 
                                 outputname='v43_can_6dists')
```

#### Panel C

Estimate a p value for the similarity between *XIST* and *KCNQ1OT1* using deduplicated GENCODE canonical lncRNAs as background and the best-fitting distribution from seekr_find_dist as the function to calculate significance.


Estimate p values for the similarity between *XIST* and all other GENCODE canonical lncRNAs.

```python
from seekr import find_dist
from seekr import find_pval
from seekr import adj_pval

v43_6dists = find_dist.find_dist(inputseq='v43_canonical.fa', k_mer=6, 
                                 log2='Log2.post', models='common10', 
                             	 subsetting=True, subset_size = 100000, 
                                 fit_model=True, statsmethod='ks',
                             	 progress_bar=True, plotfit='modelfit6.pdf', 
                                 outputname='v43_can_6dists')

xist_ot1_pval=find_pval.find_pval(seq1file='XIST.fa', 
                                  seq2file='KCNQ1OT1.fa', 
                                  mean_path='mean_6mers.npy', 
                                  std_path='std_6mers.npy',
                                  k_mer=6, fitres=v43_6dists, 
                                  log2='Log2.post', bestfit=1, 
                                  outputname='XIST_v_OT1_pval', 
                                  progress_bar=True)


v43_xist_pval=find_pval.find_pval(seq1file='v43_canonical.fa', 
                                  seq2file='XIST.fa', 
                                  mean_path='mean_6mers.npy', 
                                  std_path='std_6mers.npy',
                                  k_mer=6, fitres=v43_6dists, 
                                  log2='Log2.post', bestfit=1, 
                                  outputname='v43_v_XIST_pval', 
                                  progress_bar=True)

adjpvals=adj_pval.adj_pval(v43_xist_pval, method='fdr_bh', 
                           alpha=0.05, outputname='v43_v_XIST_pvals_bh')
```

### Figure 2 graphing with chunks

#### generate chunks

Therefore, to determine whether our select XIST-like lncRNAs harbored notable regional similarities with XIST, we fragmented the lncRNAs into ~500 nucleotide (nt) chunks and compared the chunks in each XIST-like lncRNA to each chunk within XIST.

For XIST, we fragemented it based on its sequence features: each repeat is by itself an individual chunk. For the intervals that is longer than 1000nt, we fragmented it in the same way as the other lncRNAs.

The python code for this part is available in this github repository named: generate_chunks.py
The fragmented lncRNAs and XIST fasta file are also available: v43_can500_namedchunks_500.fa and XIST_manual_chunks.fa

#### Panel B to D

Fragment select *XIST*-like lncRNAs into ~500 nucleotide (nt) chunks and compare the chunks in each *XIST*-like lncRNA to each chunk within *XIST*. Evaluate each *XIST* Repeat as a single intact chunk and separate the intervening *XIST* intervals into ~500 chunks. 

```python
from seekr.kmer_counts import BasicCounter
from seekr import find_dist
from seekr import find_pval
from seekr import kmer_heatmap
import numpy as np

v43_4dists = find_dist.find_dist(inputseq='v43_can500_namedchunks_500.fa', 
                                 k_mer=4, log2='Log2.post', models='common10',
                                 subsetting=True, subset_size = 100000, 
                                 fit_model=True, statsmethod='ks', 
                                 progress_bar=True, plotfit='modelfit4.pdf', 
                                 outputname='v43_can_4dists')


bkg_norm_4 = BasicCounter('v43_can500_namedchunks_500.fa', k=4)
bkg_norm_4.get_counts()

mean_path = 'mean_4mers.npy'
std_path = 'std_4mers.npy'
np.save(mean_path, bkg_norm_4.mean)
np.save(std_path, bkg_norm_4.std)


dlx6_xist_pval=find_pval.find_pval(seq1file='XIST_manual_chunks.fa', 
                                   seq2file='DLX6-AS1_namedchunks_500.fa',
                                   mean_path='mean_4mers.npy', 
                                   std_path='std_4mers.npy',
                                   k_mer=4, fitres=v43_4dists, 
                                   log2='Log2.post', bestfit=1, 
                                   outputname='DLX6_v_X_pvals', 
                                   progress_bar=True)

linc00632_xist_pval=find_pval.find_pval(seq1file='XIST_manual_chunks.fa', 
                                        seq2file='LINC00632_namedchunks_500.fa', 
                                        mean_path='mean_4mers.npy', 
                                        std_path='std_4mers.npy',
                                        k_mer=4, fitres=v43_4dists, 
                                        log2='Log2.post', bestfit=1,
                                        outputname='LINC00632_v_X_pvals', 
                                        progress_bar=True)

pcdh10_xist_pval=find_pval.find_pval(seq1file='XIST_manual_chunks.fa', 
                                     seq2file='PCDH10-DT_namedchunks_500.fa', 
                                     mean_path='mean_4mers.npy', 
                                     std_path='std_4mers.npy',
                                     k_mer=4, fitres=v43_4dists, 
                                     log2='Log2.post', bestfit=1, 
                                     outputname='PCDH10-DT_v_X_pvals', 
                                     progress_bar=True)


kmer_heatmap.kmer_heatmap(dlx6_xist_pval, datamin=0, datamax=1, thresh_value=0.05,
                          color_range=['#1b7837', '#ffffff', '#c51b7d'], 
                          distmetric='correlation', linkmethod='complete', 
                          hmapw_ratio=0.3, hmaph_ratio=0.3, 
                          x_tick_size=16, y_tick_size=16, 
                          cbar_font_size=16, outputname='DLX6_v_X_pvals', 
                          hformat='pdf', hdpi=300)

kmer_heatmap.kmer_heatmap(linc00632_xist_pval, datamin=0, datamax=1, thresh_value=0.05,
                          color_range=['#1b7837', '#ffffff', '#c51b7d'], 
                          distmetric='correlation', linkmethod='complete', 
                          hmapw_ratio=0.3, hmaph_ratio=0.3, 
                          x_tick_size=16, y_tick_size=16, 
                          cbar_font_size=16, outputname='LINC00632_v_X_pvals', 
                          hformat='pdf', hdpi=300)

kmer_heatmap.kmer_heatmap(pcdh10_xist_pval, datamin=0, datamax=1, thresh_value=0.05, 
                          color_range=['#1b7837', '#ffffff', '#c51b7d'], 
                          distmetric='correlation', linkmethod='complete', 
                          hmapw_ratio=0.3, hmaph_ratio=0.3, 
                          x_tick_size=16, y_tick_size=16, 
                          cbar_font_size=16, outputname='PCDH10-DT_v_X_pvals', 
                          hformat='pdf', hdpi=300)
```

#### Panel E

Plot a heatmap of *k*-mer similarity scores across fragments of *XIST*-like lncRNAs of interest. Use *k* = 3 for ease of visualizing individual *k*-mers.

```python
from seekr.kmer_counts import BasicCounter
from seekr import kmer_heatmap
import numpy as np

bkg_norm_3 = BasicCounter('v43_canonical.fa', k=6)
bkg_norm_3.get_counts()

mean_path = 'mean_3mers.npy'
std_path = 'std_3mers.npy'
np.save(mean_path, bkg_norm_3.mean)
np.save(std_path, bkg_norm_3.std)


slnc_count = BasicCounter(infasta='select_lncs.fa', 
                          outfile='select_lnc_3mers.csv', 
                          k=3, mean='mean_3mers.npy', 
                          std='std_3mers.npy', log2='Log2.post')

slnc_count. make_count_file()


kmer_heatmap.kmer_heatmap(slnc_count, datamin=0, datamax=3, thresh_value=1, 
                          color_range=['#fcfc03', '#ffffff', '#1907e6'], 
                          distmetric='correlation', linkmethod='complete', 
                          hmapw_ratio=0.3, hmaph_ratio=0.3, 
                          x_tick_size=16, y_tick_size=16, 
                          cbar_font_size=16, outputname='select_lnc_3mers', 
                          hformat='pdf', hdpi=300)
```

#### Panel F

Plot a barplot of most-enriched *k*-mers across fragments of *XIST*-like lncRNAs of interest. Use *k* = 3 for ease of visualizing individual *k*-mers.

```python
from seekr.kmer_counts import BasicCounter
from seekr import kmer_msd_barplot
import numpy as np

bkg_norm_3 = BasicCounter('v43_canonical.fa', k=6)
bkg_norm_3.get_counts()

mean_path = 'mean_3mers.npy'
std_path = 'std_3mers.npy'
np.save(mean_path, bkg_norm_3.mean)
np.save(std_path, bkg_norm_3.std)


kmer_msd_barplot.kmer_msd_barplot(inputfile='select_lncs.fa', 
                                  mean='mean_3mers.npy', std='std_3mers.npy', 
                                  log2 = 'Log2.post', k=3, 
                                  sortstat='mean', sortmethod='descending', 
                                  topkmernumber=10, xlabelsize=20, ylabelsize=20, 
                                  xticksize=20, yticksize=20, 
                                  outputname='select_lncs_msd_barplot3', 
                                  pformat='pdf', pdpi=300)
```

#### Panel G

Plot a barplot of most variable *k*-mers across fragments of *XIST*-like lncRNAs of interest. Use *k* = 3 for ease of visualizing individual *k*-mers.

```python
from seekr.kmer_counts import BasicCounter
from seekr import kmer_count_barplot
import numpy as np

bkg_norm_3 = BasicCounter('v43_canonical.fa', k=6)
bkg_norm_3.get_counts()

mean_path = 'mean_3mers.npy'
std_path = 'std_3mers.npy'
np.save(mean_path, bkg_norm_3.mean)
np.save(std_path, bkg_norm_3.std)

kmer_count_barplot.kmer_count_barplot(inputfile='select_lncs.fa', 
                                      mean='mean_3mers.npy', std='std_3mers.npy',
                                      log2='Log2.post', k=3, 
                                      sortmethod='descending', topkmernumber=10,
                                      xlabelsize=20, ylabelsize=20, 
                                      xticksize=20, yticksize=20, legendsize=12,
                                      outputname='select_lncs_barplot3', 
                                      pformat='pdf', pdpi=300)
```

#### Panel H

Compare *k*-mer similarity among fragments of *XIST* at *k*-mer length *k* = 4. Calculate p-values and plot them in a heatmap.

```python
from seekr.kmer_counts import BasicCounter
from seekr import find_pval
from seekr import kmer_heatmap
import numpy as np

bkg_norm_4 = BasicCounter('v43_can500_namedchunks_500.fa', k=4)
bkg_norm_4.get_counts()

mean_path = 'mean_4mers.npy'
std_path = 'std_4mers.npy'
np.save(mean_path, bkg_norm_4.mean)
np.save(std_path, bkg_norm_4.std)


xist_xist_pval=find_pval.find_pval(seq1file='XIST_manual_chunks.fa', 
                                   seq2file='XIST_manual_chunks.fa', 
                                   mean_path='mean_4mers.npy', 
                                   std_path='std_4mers.npy',
                                   k_mer=4, fitres=v43_4dists, 
                                   log2='Log2.post', bestfit=1, 
                                   outputname='X_v_X_pvals',
                                   progress_bar=True)


kmer_heatmap.kmer_heatmap(xist_xist_pval, datamin=0, datamax=1, thresh_value=0.05,
                          color_range=['#1b7837', '#ffffff', '#c51b7d'], 
                          distmetric='correlation', linkmethod='complete', 
                          hmapw_ratio=0.4, hmaph_ratio=0.3, 
                          x_tick_size=16, y_tick_size=16, 
                          cbar_font_size=16, outputname='X_v_X_pvals', 
                          hformat='pdf', hdpi=300)
```

#### Panel I

Compare *k*-mer similarity among fragments of *XIST* at *k*-mer length *k* = 4. From p-values calculated above, create a dendrogram of fragment similarities across *XIST*.

```python
from seekr import kmer_dendrogram

kmer_dendrogram.kmer_dendrogram(xist_xist_pval, dendro_direct='column', 
                                distmetric='correlation', linkmethod='complete', 
                                plot_ht=8, wd_ratio=0.5, leaf_font_size = 16, 
                                outputname='X_v_X_dendrogram', 
                                pformat='pdf', pdpi=300)
```

#### Panel S1

Compare *k*-mer similarity among fragments of *XIST* at *k*-mer length *k* = 4 and use Leiden to assign fragments to communities.

```python
from seekr.kmer_counts import BasicCounter
from seekr import kmer_leiden
import numpy as np

bkg_norm_4 = BasicCounter('v43_can500_namedchunks_500.fa', k=4)
bkg_norm_4.get_counts()

mean_path = 'mean_4mers.npy'
std_path = 'std_4mers.npy'
np.save(mean_path, bkg_norm_4.mean)
np.save(std_path, bkg_norm_4.std)


kmer_leiden.kmer_leiden(inputfile='XIST_manual_chunks.fa', 
                        mean='mean_4mers.npy', std='std_4mers.npy', 
                        k=4, algo='RBERVertexPartition', 
                        rs=1.1, pearsoncutoff=0.15, 
                        csvfile='XIST_manual_chunks')
```

The result files: XIST_manual_chunks_nodes_leiden.csv and XIST_manual_chunks_edges_leiden.csv are then directly imported into Gephi as nodes and edges files for network plotting with Force Atalas2 layout.

### Table 2 and Table S2

Fragment GENCODE canonical lncRNAs into ~500nt chunks and compare each GENCODE lncRNA chunk to each *XIST* Repeat, using *k*-mer *k* = 4. Calculate p values of similarity using GENCODE lncRNA chunks as background and fitting to a lognormal distribution.

Define similarity as any fragment whose SEEKR-derived similarity to an *XIST* Repeat has a p value of less than 0.05. Tally the number of *XIST*-similar fragments in each lncRNA and used these sums to rank lncRNAs by their overall *XIST*-likeness.

XIST_repeats.fa is also provided under this github repository.

```python
from seekr.kmer_counts import BasicCounter
from seekr import find_pval
import numpy as np

bkg_norm_4 = BasicCounter('v43_can500_namedchunks_500.fa', k=4)
bkg_norm_4.get_counts()

mean_path = 'mean_4mers.npy'
std_path = 'std_4mers.npy'
np.save(mean_path, bkg_norm_4.mean)
np.save(std_path, bkg_norm_4.std)


v43_ABDEF_pval=find_pval.find_pval(seq1file='v43_can500_namedchunks_500.fa', 
                                   seq2file='XIST_repeats.fa', 
                                   mean_path='mean_4mers.npy', 
                                   std_path='std_4mers.npy',
                                   k_mer=4, fitres=v43_4dists, 
                                   log2='Log2.post', bestfit=1, 
                                   outputname='v43_chunks_v_ABDEF_pvals', 
                                   progress_bar=True)
```

The result file v43_chunks_v_ABDEF_pvals.csv is then processed in R with the code: chunkABDE_organize.R to produce the count file: v43_lncRNA_ABDEF_count.csv. This file is then used to generate Table 2 and Table S2


### coPARSE analysis

All codes for this part of analysis is saved under the **coparse folder**.

From [GENCODE](https://www.gencodegenes.org/) download lncRNA fasta and gtf file for human: v25 and mouse: vM10.

Filter both the human and mouse sequences and setting the background distribution of seekr r values using filtered human lncRNAs. Exract human and mouse sequences corresponds to Table S4 and Table S5a in the [coPARSE paper](https://www.nature.com/articles/s41588-023-01620-7). Calculate the seekr r value, p value and adjusted p value for all human vs mouse pair in Table S4 and S5a. Codes for this part are saved as **coparse_analysis.py**, which can be downloaded from this repo. Here we did not use the command line function as there are functions used that are outside the seekr package and are not command line compatible. So it is easier to accomplish the whole session in Python. 

Plot seekr r value background distribution as a density plot; overlap with the human vs mouse r values in Table S4 and Table S5a as histogram. The histogram is then colored according to their p values. Codes for this part are saved as **coparse_plotting.R**. 

Cut human and mouse sequences in Table S4 and Table S5a into chunks, using code: **coparse_generate_chunks.py**.

Calculate the seekr r value, p value and adjusted p value for sequences in Table S4 and S5a vs XIST full sequence; and chunked sequences of Table S4 and Table S5 vs XIST repeats. All functions needed here are within the seekr package. So here we provide the codes as in command line: **coparse_XIST_S45.sh**.

Analyze the results to get the significance chunk count per XIST repeat per gene and then combine all results together. Codes are saved as **coparse_chunk_analysis.R**.


### randomness analysis

All codes for this part of analysis is saved under the **randomness folder**.

Using previously generated data, run v43_canonical.fa against XIST.fa at k=6 with different subset size (10k, 100k, 1M and all-no subsetting) and fit the 'common10' models. Each simulate 10 times. Codes used here are saved as **v43_XIST_fit_randomness_10k/100k/1M/full.sh**.

From the fitted model files, get the best fit model name and calculate the significant gene counts by either p value <0.05 or adjusted p value<0.05 corresponds to each simulation. Plot the significante gene counts against different subset size using codes: **randomness_boxplot.R**.

Use subset size = 100k, fit the data with 'all' models. Calcualte the p values and adjusted p values for all model fit, using the best fitted model, the second best fitted model and the third best fitted model. Codes saved as: **v43_XIST_fit_randomness_100k_all.sh**.

Compare this results with the previous 'common10' best fitted model. Generate the dataframe that includes the r value that corresponds to 95% percentile of each fitted model (where p=0.05), and the significant gene count using p value<0.05 and adjusted p value<0.05. Codes saved as **randomness_allvscommon10.py**. And plot the boxplot of the comparison results using **randomness_boxplot.R**.


### human (s)eCLIP analysis

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

