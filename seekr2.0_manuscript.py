import numpy as np
from seekr.kmer_counts import BasicCounter
from seekr.pearson import pearson as seekrPearson
from seekr import fasta
from seekr import filter_gencode
from seekr import find_dist
from seekr import find_pval
from seekr import adj_pval
from seekr import kmer_heatmap
from seekr import kmer_msd_barplot
from seekr import kmer_count_barplot
from seekr import kmer_dendrogram
from seekr import kmer_leiden



downloader = fasta.Downloader()
downloader.get_gencode(biotype='lncRNA', species='human', gtf=True, release='43', unzip=True)

headers, seqs = filter_gencode.filter_gencode(fasta_path='v43_lncRNA.fa', 
                                              gtf_path='v43_lncRNA.chr_patch_hapl_scaff.annotation.gtf',
                                              len_threshold=500, canonical=True,
                                              rm_dup=True, outputname='v43_canonical')


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




kmer_msd_barplot.kmer_msd_barplot(inputfile='select_lncs.fa', 
                                  mean='mean_3mers.npy', std='std_3mers.npy', 
                                  log2 = 'Log2.post', k=3, 
                                  sortstat='mean', sortmethod='descending', 
                                  topkmernumber=10, xlabelsize=20, ylabelsize=20, 
                                  xticksize=20, yticksize=20, 
                                  outputname='select_lncs_msd_barplot3', 
                                  pformat='pdf', pdpi=300)


kmer_count_barplot.kmer_count_barplot(inputfile='select_lncs.fa', 
                                      mean='mean_3mers.npy', std='std_3mers.npy',
                                      log2='Log2.post', k=3, 
                                      sortmethod='descending', topkmernumber=10,
                                      xlabelsize=20, ylabelsize=20, 
                                      xticksize=20, yticksize=20, legendsize=12,
                                      outputname='select_lncs_barplot3', 
                                      pformat='pdf', pdpi=300)



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



kmer_dendrogram.kmer_dendrogram(xist_xist_pval, dendro_direct='column', 
                                distmetric='correlation', linkmethod='complete', 
                                plot_ht=8, wd_ratio=0.5, leaf_font_size = 16, 
                                outputname='X_v_X_dendrogram', 
                                pformat='pdf', pdpi=300)



kmer_leiden.kmer_leiden(inputfile='XIST_manual_chunks.fa', 
                        mean='mean_4mers.npy', std='std_4mers.npy', 
                        k=4, algo='RBERVertexPartition', 
                        rs=1.1, pearsoncutoff=0.15, 
                        csvfile='XIST_manual_chunks')



v43_ABDEF_pval=find_pval.find_pval(seq1file='v43_can500_namedchunks_500.fa', 
                                   seq2file='XIST_repeats.fa', 
                                   mean_path='mean_4mers.npy', 
                                   std_path='std_4mers.npy',
                                   k_mer=4, fitres=v43_4dists, 
                                   log2='Log2.post', bestfit=1, 
                                   outputname='v43_chunks_v_ABDEF_pvals', 
                                   progress_bar=True)