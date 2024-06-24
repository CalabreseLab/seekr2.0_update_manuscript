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


# cut seqs into chunks
v25s4 = 'v25_S4.fa'
v25s4_div_seqs, v25s4_div_seq_headers = divide_fasta_to_chunks(v25s4, 500, 1)

# save these gencode chunks into a fasta file
save_seqs_to_fasta(v25s4_div_seqs, v25s4_div_seq_headers, 'v25_S4_500chunk.fa')


# cut seqs into chunks
v25s5 = 'v25_S5.fa'
v25s5_div_seqs, v25s5_div_seq_headers = divide_fasta_to_chunks(v25s5, 500, 1)

# save into a fasta file
save_seqs_to_fasta(v25s5_div_seqs, v25s5_div_seq_headers, 'v25_S5_500chunk.fa')


vm10s4 = 'vM10_S4.fa'
vm10s4_div_seqs, vm10s4_div_seq_headers = divide_fasta_to_chunks(vm10s4, 500, 1)

# save  into a fasta file
save_seqs_to_fasta(vm10s4_div_seqs, vm10s4_div_seq_headers, 'vM10_S4_500chunk.fa')

vm10s5 = 'vM10_S5.fa'
vm10s5_div_seqs, vm10s5_div_seq_headers = divide_fasta_to_chunks(vm10s5, 500, 1)

# save  into a fasta file
save_seqs_to_fasta(vm10s5_div_seqs, vm10s5_div_seq_headers, 'vM10_S5_500chunk.fa')



