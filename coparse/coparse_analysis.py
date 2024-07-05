
import numpy as np
import pandas as pd


from seekr.kmer_counts import BasicCounter 
from seekr.pearson import pearson as seekrPearson
from seekr.fasta_reader import Reader as seekrReader

from seekr.find_dist import find_dist
from seekr.find_pval import find_pval
from seekr.adj_pval import adj_pval

from seekr.filter_gencode import filter_gencode

#####################
## extract sequences based on S4 and S5 tables
#####################

# filter the gencode fasta file to only keep seqs >=500nt and remove duplicates

headers, seqs = filter_gencode(fasta_path='gencode.v25.lncRNA_transcripts.fa', 
                               gtf_path='gencode.v25.long_noncoding_RNAs.gtf',
                               len_threshold=500, canonical=False,
                               rm_dup=True, outputname='v25_500_rd')

headers, seqs = filter_gencode(fasta_path='gencode.vM10.lncRNA_transcripts.fa', 
                               gtf_path='gencode.vM10.long_noncoding_RNAs.gtf',
                               len_threshold=500, canonical=False,
                               rm_dup=True, outputname='vM10_500_rd')

# fit filtered human vs filter human and generate the background model
v25_6dists = find_dist(inputseq='v25_500_rd.fa', k_mer=6, 
                       log2='Log2.post', models='common10', 
                       subsetting=True, subset_size = 100000, 
                       fit_model=True, statsmethod='ks',
                       progress_bar=True, plotfit='modelfit6.pdf',
                       outputname='v25_can_6dist')


# extract seqs based on S4 table
from Bio import SeqIO

# define function to extract seqs based on transcript ids
def extract_seqs(fasta_file, output_fasta, transcript_ids):
    # Read all sequences from the FASTA file into a list of tuples for comprehensive access
    sequences = [(seq.id.split('|')[0].split('.')[0], seq) for seq in SeqIO.parse(fasta_file, 'fasta')]

    # Create an output list for storing sequences in the order they appear in the CSV
    output_seqs = []

    # Append sequences to the output list based on their appearance in the transcript_ids list
    for tid in transcript_ids:
        matches = [seq for key, seq in sequences if key == tid]
        if len(matches) == 0:
            print(f"No matches found for TID: {tid}")
        elif len(matches) > 1:
            print(f"Multiple matches found for TID: {tid}")
        
        output_seqs.extend(matches)

    # Write the ordered sequences to a new FASTA file, including duplicates, with original headers
    SeqIO.write(output_seqs, output_fasta, 'fasta')


# extract human S4 seqs
# Load CSV file
csv_file = 'coparse_tableS4.csv'
df = pd.read_csv(csv_file)

# Load all transcript IDs to maintain order and duplicates
h_transcript_ids = df['human transcript'].tolist()

# Load FASTA file
fasta_file = 'gencode.v25.lncRNA_transcripts.fa'
output_fasta = 'v25_S4.fa'

# extract seqs
extract_seqs(fasta_file, output_fasta, h_transcript_ids)


# extract mouse S4 seqs
# Load all transcript IDs to maintain order and duplicates
m_transcript_ids = df['mouse transcript'].tolist()

# Load FASTA file
fasta_file = 'gencode.vM10.lncRNA_transcripts.fa'
output_fasta = 'vM10_S4.fa'

# extract seqs
extract_seqs(fasta_file, output_fasta, m_transcript_ids)


# extract human S5 seqs
# Load CSV file
df = pd.read_csv('coparse_tableS5a.csv')

# Load all transcript IDs to maintain order and duplicates
h_transcript_ids = df['human transcript'].tolist()
# Load FASTA file
fasta_file = 'gencode.v25.lncRNA_transcripts.fa'
output_fasta = 'v25_S5a.fa'

# extract seqs
extract_seqs(fasta_file, output_fasta, h_transcript_ids)

# extract mouse S5 seqs
# Load all transcript IDs to maintain order and duplicates
m_transcript_ids = df['mouse transcript'].tolist()
# Load FASTA file
fasta_file = 'gencode.vM10.lncRNA_transcripts.fa'
#fasta_file = 'v25_500_rd.fa'
output_fasta = 'vM10_S5a.fasta'

# extract seqs
extract_seqs(fasta_file, output_fasta, m_transcript_ids)

#######################
#get the seerk r value for human vs human and human vs mouse in S4 and S5
###########################

# set the normalization vectors for seekr
humanbkg = BasicCounter('v25_500_rd.fa', k=6)
humanbkg.get_counts()
mean_path = 'v25_mean6.npy'
std_path = 'v25_std6.npy'
np.save(mean_path, humanbkg.mean)
np.save(std_path, humanbkg.std)

# Count k-mers of human lnc
human_counter = BasicCounter('v25_500_rd.fa',outfile='human_6mer.npy',mean='v25_mean6.npy',std='v25_std6.npy',k=6)
human_counter.make_count_file()

# Find similarities of human vs human
hsim = seekrPearson(human_counter.counts,human_counter.counts,outfile='human_vs_human_6mer_fullrvals.npy')

# get the values in upper triangle of gc_sim and saved as 1D np array
hsim_triu = hsim[np.triu_indices(hsim.shape[0], k=1)]

# save gc_sim_triu as feature file
hsim_triu_df=pd.DataFrame(hsim_triu, columns=['Data'])
# hsim_triu_df.to_feather('human_vs_human_6mer_rvals.feather')

# randomly sample 1M data points and save as csv
hh_ss=hsim_triu_df.sample(n=1000000, replace=False, random_state=42)
hh_ss.to_csv('human_vs_human_6mer_rvals_1m.csv',index=False)
# this is used as the background distribution density plot in the figure

# Count k-mers of S4 human and mouse
s4h_counter = BasicCounter('v25_S4.fasta',outfile='human_s4_6mer.npy',mean='v25_mean6.npy',std='v25_std6.npy',k=6)
s4m_counter = BasicCounter('vM10_S4.fasta',outfile='mouse_s4_6mer.npy',mean='v25_mean6.npy',std='v25_std6.npy',k=6)

s4h_counter.make_count_file()
s4m_counter.make_count_file()

# Find similarities
s4sim = seekrPearson(s4h_counter.counts,s4m_counter.counts,outfile='s4_fullrvals.npy')

# only use the diagonal which is the r values of the human vs mouse pairs in table S4
s4sim_dia=np.diag(s4sim)
s4sim_df=pd.DataFrame(s4sim_dia)
s4sim_df.to_csv('s4_diag_rvals.csv',index=False)
# this is the histogram in the figure

# Count k-mers of S5 human and mouse
s5h_counter = BasicCounter('v25_S5a.fasta',outfile='human_s4_6mer.npy',mean='v25_mean6.npy',std='v25_std6.npy',k=6)
s5m_counter = BasicCounter('vM10_S5a.fasta',outfile='mouse_s4_6mer.npy',mean='v25_mean6.npy',std='v25_std6.npy',k=6)

s5h_counter.make_count_file()
s5m_counter.make_count_file()

# Find similarities
s5sim = seekrPearson(s5h_counter.counts,s5m_counter.counts,outfile='s5_fullrvals.npy')

# only use the diagonal which is the r values of the human vs mouse pairs in table S5
s5sim_dia=np.diag(s5sim)
s5sim_df=pd.DataFrame(s5sim_dia)
s5sim_df.to_csv('s5_diag_rvals.csv',index=False)
# this is the histogram in the figure

# check sig human mouse pair compared to human background dist in S4
s4pvals=find_pval(seq1file='v25_S4.fasta', seq2file='vM10_S4.fasta', mean_path='v25_mean6.npy', std_path='v25_std6.npy',
                  k_mer=6, fitres=v25_6dists, log2='Log2.post',
                  bestfit=1, outputname='s4_pvals', progress_bar=True)

s4adjpvals=adj_pval(s4pvals, method='fdr_bh', alpha=0.05, outputname='s4_bh_adjpval')

# check sig human mouse pair compared to human background dist in S5
s5pvals=find_pval(seq1file='v25_S5a.fasta', seq2file='vM10_S5a.fasta', mean_path='v25_mean6.npy', std_path='v25_std6.npy',
                  k_mer=6, fitres=v25_6dists, log2='Log2.post',
                  bestfit=1, outputname='s5a_pvals', progress_bar=True)

s5adjpvals=adj_pval(s5pvals, method='fdr_bh', alpha=0.05, outputname='s5_bh_adjpval')

# only get the diagonal values as these are the p vals and adjusted p vals for the human vs mouse pairs in S4 and S5
s4pvals_dia=np.diag(s4pvals)
s4pvals_df=pd.DataFrame(s4pvals_dia)
s4pvals_df.to_csv('s4_diag_pvals.csv',index=False)

s4adjpvals_dia=np.diag(s4adjpvals)
s4adjpvals_df=pd.DataFrame(s4adjpvals_dia)
s4adjpvals_df.to_csv('s4_diag_bh_adjpvals.csv',index=False)

s5pvals_dia=np.diag(s5pvals)
s5pvals_df=pd.DataFrame(s5pvals_dia)
s5pvals_df.to_csv('s5_diag_pvals.csv',index=False)

s5adjpvals_dia=np.diag(s5adjpvals)
s5adjpvals_df=pd.DataFrame(s5adjpvals_dia)
s5adjpvals_df.to_csv('s5_diag_bh_adjpvals.csv',index=False)
# these are used for coloring of the histograms


# Count k-mers of manually selected highly conserved human and mouse seqs 
human_hc_counter = BasicCounter('human_highly_conserved.fa',outfile='human_hc_6mer.npy',mean='v25_mean6.npy',std='v25_std6.npy',k=6)
mouse_hc_counter = BasicCounter('mouse_highly_conserved.fa',outfile='mouse_hc_6mer.npy',mean='v25_mean6.npy',std='v25_std6.npy',k=6)

human_hc_counter.make_count_file()
mouse_hc_counter.make_count_file()

# Find similarities
hcsim = seekrPearson(human_hc_counter.counts,mouse_hc_counter.counts,outfile='hc_hvsm.npy')

hcsim_dia=np.diag(hcsim)
hcsim_df=pd.DataFrame(hcsim_dia)
hcsim_df.to_csv('highly_conserved_rvals.csv',index=False)
# these are the vertical lines annotated in the figure

# plotting codes are in coparse_plotting.R









