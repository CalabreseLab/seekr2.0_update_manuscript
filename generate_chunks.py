# generate chunks for input fasta
# save gene name and corrdinates for each chunk in the header

import numpy as np

from seekr.fasta_reader import Reader as seekrReader


# divide sequence into chunks around a target length without overlap
# length of seq divided by target length, the quotient is the number of chunks
# if remainder is greater than flex_length, add one more chunk
# if remainder is less than flex_length, the total num of chunk is the quotient
# divide the sequence equally into the total num of chunks, the remiander goes to the last chunk
def divide_seq(seq, target_length, flex_length): 
    seq_length = len(seq)
    quotient = seq_length // target_length
    remainder = seq_length % target_length
    if remainder > flex_length:
        num_chunk = quotient + 1
    else:
        num_chunk = quotient
    chunk_length = seq_length // num_chunk
    chunks = [seq[i:i+chunk_length] for i in range(0, seq_length, chunk_length)]
    # save the index range for each chunk
    # if last index is greater than seq_length, save seq_length instead
    chunk_index = [(i, np.min([i+chunk_length, seq_length])) for i in range(0, seq_length, chunk_length)]
    
    # append the last chunk to the second to the last chunk
    # if the last chunk length is less than chunk_length
    if len(chunks[-1]) < chunk_length:
        chunks[-2] = chunks[-2] + chunks[-1]
        # update chunk_index accordingly
        chunk_index[-2] = (chunk_index[-2][0], chunk_index[-1][1])
        del chunks[-1]
        del chunk_index[-1]
    return chunks, chunk_index


# read in fasta file and divide all seqs using divide_seq function
# write each divided chunks into a list

def divide_fasta_to_chunks(seqpath, target_length, flex_length): 
    seqs = seekrReader(seqpath).get_seqs()
    headers=seekrReader(seqpath).get_headers()
    div_seqs = []
    div_seq_headers = []
    for i, seq in enumerate(seqs):
        # only keep seqs longer than flex_length
        if len(seq) > flex_length:
            # if seq is longer than target_length, divide it into chunks
            if len(seq) > target_length:
                chunks, chunk_index = divide_seq(seq, target_length, flex_length)
                div_seqs.append(chunks)
                # save the header for all chunks
                chunk_headers=[headers[i] + '_' + str(index[0]) + '_' + str(index[1]) for index in chunk_index]
                div_seq_headers.append(chunk_headers)
            # if seq is shorter than target_length, append it to the list
            else:
                div_seqs.append([seq])
                # add seq length to the header
                div_seq_headers.append([headers[i] + '_0_' + str(len(seq))])
    # flatten the list of lists into a list
    div_seqs = [item for sublist in div_seqs for item in sublist]
    div_seq_headers = [item for sublist in div_seq_headers for item in sublist]
    return div_seqs, div_seq_headers


# save all element of a list of seqs into one fasta file
# close the files after written to avoid too many open files error

def save_seqs_to_fasta(seqs, seq_names, seqpath):
    seqfile = open(seqpath, 'w')
    for i, seq in enumerate(seqs):
        # save fasta name from corresponding element in seq_names
        seqfile.write(seq_names[i] + '\n' + seq + '\n')
    seqfile.close()


# cut canonical v43 seqs into chunks
gencode = 'v43_canonical.fa'
gencode_div_seqs, gencode_div_seq_headers = divide_fasta_to_chunks(gencode, 500, 200)

# save these gencode chunks into a fasta file
save_seqs_to_fasta(gencode_div_seqs, gencode_div_seq_headers, 'v43_can500_namedchunks_500.fa')


# divide DLX6, LINC and PCDH10 into chunks
dlx6 = 'ENST00000430027.3_(DLX6-AS1).fa'
dlx6_div_seqs, dlx6_div_seq_headers = divide_fasta_to_chunks(dlx6, 500, 200)
# 42150 chunks in total before

# save  into a fasta file
save_seqs_to_fasta(dlx6_div_seqs, dlx6_div_seq_headers, 'DLX6-AS1_namedchunks_500.fa')

linc = 'ENST00000648200.2_(LINC00632).fa'
linc_div_seqs, linc_div_seq_headers = divide_fasta_to_chunks(linc, 500, 200)
# 42150 chunks in total before

# save  into a fasta file
save_seqs_to_fasta(linc_div_seqs, linc_div_seq_headers, 'LINC00632_namedchunks_500.fa')

pcdh10 = 'ENST00000667505.1_(PCDH10-DT).fa'
pcdh10_div_seqs, pcdh10_div_seq_headers = divide_fasta_to_chunks(pcdh10, 500, 200)
# 42150 chunks in total before

# save  into a fasta file
save_seqs_to_fasta(pcdh10_div_seqs, pcdh10_div_seq_headers, 'PCDH10-DT_namedchunks_500.fa')



