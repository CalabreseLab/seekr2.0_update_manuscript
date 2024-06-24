# generate bedfiles for sig chunks by converting to genomic coordinates
# randomize the chunks and generate the bedfiles in the same way

# set working directory to where the files are located

library(Biostrings) # needed to be loaded first, before rtracklayer
library(rtracklayer)

###########################################
# Function to map transcript coordinates to genomic coordinates
# for pos strand
map_transcript_to_genomic_p <- function(transcript_start, transcript_end, texons) {
  # Find the exons that the transcript coordinates overlap
  start_exon = which(transcript_start <= texons$cumulative_length)[1]
  end_exon = which(transcript_end <= texons$cumulative_length)[1]
  
  # Initialize variables to keep track of genomic coordinates
  genomic_coordinates <- list()
  
  # Loop over exons from start_exon to end_exon
  for (i in start_exon:end_exon) {
    if (i == start_exon) {
      # Calculate the genomic start position for the first overlapping exon
      genomic_start = texons$start[i] + (transcript_start - (ifelse(i == 1, 0, texons$cumulative_length[i - 1])) - 1)
    } else {
      genomic_start = texons$start[i]
    }
    
    if (i == end_exon) {
      # Calculate the genomic end position for the last overlapping exon
      genomic_end = texons$start[i] + (transcript_end - (ifelse(i == 1, 0, texons$cumulative_length[i - 1])) - 1)
    } else {
      genomic_end = texons$end[i]
    }
    
    # Append the genomic coordinate segment
    genomic_coordinates[[length(genomic_coordinates) + 1]] <- c(genomic_start, genomic_end)
  }
  
  # Convert list to matrix for easier viewing
  do.call(rbind, genomic_coordinates)
}

# Function to map transcript coordinates to genomic coordinates
# for neg strand
map_transcript_to_genomic_n <- function(transcript_start, transcript_end, texons) {
  # Find the exons that the transcript coordinates overlap
  start_exon = which(transcript_start <= texons$cumulative_length)[1]
  end_exon = which(transcript_end <= texons$cumulative_length)[1]
  
  # Initialize variables to keep track of genomic coordinates
  genomic_coordinates <- list()
  
  # Loop over exons from start_exon to end_exon
  for (i in start_exon:end_exon) {
    if (i == start_exon) {
      # Calculate the genomic start position for the first overlapping exon
      # going in reverse 
      genomic_start = texons$end[i] - (transcript_start - (ifelse(i == 1, 0, texons$cumulative_length[i - 1])) - 1)
    } else {
      genomic_start = texons$end[i]
    }
    
    if (i == end_exon) {
      # Calculate the genomic end position for the last overlapping exon
      genomic_end = texons$end[i] - (transcript_end - (ifelse(i == 1, 0, texons$cumulative_length[i - 1])) - 1)
    } else {
      genomic_end = texons$start[i]
    }
    
    # Append the genomic coordinate segment
    genomic_coordinates[[length(genomic_coordinates) + 1]] <- c(genomic_start, genomic_end)
  }
  
  # Convert list to matrix for easier viewing
  do.call(rbind, genomic_coordinates)
}

# function to take in the chunks and exons info and convert to bedfile
convert_transcript_to_genome <- function(sigc, exons) {
  
  sigbed<-data.frame(chrom=character(),
                     chromStart=numeric(),
                     chromEnd=numeric(),
                     name=character(),
                     score=numeric(),
                     strand=character())
  
  
  for (n in 1:nrow(sigc)) {
    
    temp_chunk<-sigc[n,]
    
    temp_gene<-exons[which(exons$transcript_id==temp_chunk$tid & exons$gene_id==temp_chunk$gid),]
    
    if (nrow(temp_gene)==0) {
      
      print('no matched info in gtf')
      print(paste0('transcriptID:',temp_chunk$tid))
      print(paste0('name:',temp_chunk$name))
      break
      
    } else {
      temp_chunk_chr<-as.character(unique(temp_gene$seqnames))
      temp_chunk_name<-paste0(temp_chunk$tid,'_',temp_chunk$tstart)
      temp_chunk_strand<-as.character(unique(temp_gene$strand))
      # convert python 0 start to 1 start
      temp_tstart<-as.numeric(temp_chunk$tstart)+1
      temp_tend<-as.numeric(temp_chunk$tend)
      # the order is right directly from the dataframe
      # Calculate the cumulative transcript length at the end of each exon
      temp_gene$cumulative_length = cumsum(temp_gene$width)
      
      if (unique(temp_gene$strand)=='+') {
        
        genomic_coords <- map_transcript_to_genomic_p(temp_tstart, temp_tend, temp_gene)
        genomic_coords<-as.data.frame(genomic_coords)
        colnames(genomic_coords)<-c('chromStart','chromEnd')
        genomic_coords$chrom<-temp_chunk_chr
        genomic_coords$name<-temp_chunk_name
        genomic_coords$score<-0
        genomic_coords$strand<-temp_chunk_strand
        
        
        sigbed<-rbind(sigbed,genomic_coords[,c('chrom','chromStart','chromEnd','name','score','strand')])
        
        
      } else if (unique(temp_gene$strand)=='-') {
        
        genomic_coords <- map_transcript_to_genomic_n(temp_tstart, temp_tend, temp_gene)
        genomic_coords<-as.data.frame(genomic_coords)
        # reverse the coords back to normal start < end
        colnames(genomic_coords)<-c('chromEnd','chromStart')
        genomic_coords$chrom<-temp_chunk_chr
        genomic_coords$name<-temp_chunk_name
        genomic_coords$score<-0
        genomic_coords$strand<-temp_chunk_strand
        
        
        sigbed<-rbind(sigbed,genomic_coords[,c('chrom','chromStart','chromEnd','name','score','strand')])
        
        
      } else {
        print('cannot decide which strand the gene belongs to')
        print(paste0('transcriptID:',temp_chunk$tid))
        print(paste0('name:',temp_chunk$name))
        break
      }
      
    }
    
  }
  
  return(sigbed)
  
}

##################################################################
# import table with k_tot_1 and k_tot_2
exp<-read.csv('Table_S2_v43_lncRNA_ABDEF_count_expression_4_23_2024_values.csv',header=T)
colnames(exp)[1]<-'geneID'

exp$k_tot_mean<-(exp$K_tot_1+exp$K_tot_2)/2
# filter for k_tot_mean >0.0625

exp.fl<-exp[which(exp$k_tot_mean>0.0625),]

# this is the list of genes that we can accept
genelist<-exp.fl$geneID
genelist<-paste0(genelist,'|')

# extract info from sig chunks
chunks<-read.csv('v43_chunks_v_ABDEF_pvals.csv', header=T)

colnames(chunks)[1]<-'chunkID'

# filter the chunks list for chunks within the genelist
chunknames<-strsplit(chunks$chunkID,'_',fixed=T)
sum(lengths(chunknames)!=3) # check length
chunks$geneID<-sapply(chunknames,'[[',1)

expchunks<-chunks[(chunks$geneID %in% genelist),]

###########
# get gtf info
gtf<-import('gencode.v43.long_noncoding_RNAs.gtf')

# get the transcript id, gene id, strand, start and end coordinates for each exon
exons <- gtf[ mcols(gtf)$type == "exon" ]

exons<-as.data.frame(exons)

##################################
# filter for sig chunks
# run this code for each repeat separately and then run all the following code to generate

rpvec<-c('rA','rF','rB1','rB2','rD','rE')
expvec<-c('rA_mc','rFbroad','rB1','rB2','rD','rE_mc')

# loop through each repeat to create the sig chunk bedfile and its corresponding randomized bedfile

for (m in 1:length(rpvec)) {

  rp<-rpvec[m]
  exp<-expvec[m]
  sigchk<-expchunks$chunkID[which(expchunks[[exp]]<0.05)]


  #####################################

  sigchk<-strsplit(sigchk,'|',fixed=T)
  sum(lengths(sigchk)!=8) # check length

  sigc_tid<-sapply(sigchk,'[[',1)
  sigc_gid<-sapply(sigchk,'[[',2)
  sigc_name<-sapply(sigchk,'[[',6)
  sigc_tlen<-as.numeric(sapply(sigchk,'[[',7))
  sigc_coord<-sapply(sigchk,'[[',8)

  sigc_coord<-strsplit(sigc_coord,'_',fixed=T)
  sigc_tstart<-as.numeric(sapply(sigc_coord,'[[',2))
  sigc_tend<-as.numeric(sapply(sigc_coord,'[[',3))

  sigc<-as.data.frame(cbind(sigc_tid,sigc_gid,sigc_name))
  colnames(sigc)<-c('tid','gid','name')
  sigc$tlen<-sigc_tlen
  sigc$tstart<-sigc_tstart
  sigc$tend<-sigc_tend


  #############################################
  # extract info from gtf

  chunk_bedfile<-convert_transcript_to_genome(sigc, exons)

  write.table(chunk_bedfile,paste0('xist_sigchunk_',rp,'_1.bed'),sep='\t', quote=F, col.names = F, row.names = F)

  ###########################################################
  # randomize the chunks among the transcripts level

  # parse the headers to get tid and gid and name
  # construct something similar to the sigc
  pool_headers<-strsplit(genelist,'|',fixed=T)
  sum(lengths(pool_headers)!=7)

  pool_tid<-sapply(pool_headers,'[[',1)
  pool_gid<-sapply(pool_headers,'[[',2)
  pool_name<-sapply(pool_headers,'[[',6)
  pool_len<-as.numeric(sapply(pool_headers,'[[',7))

  pool_comb<-as.data.frame(cbind(pool_tid,pool_gid,pool_name))
  colnames(pool_comb)<-c('tid','gid','name')
  pool_comb$tlen<-pool_len
  pool_comb$tstart<-0
  # this is the pool to randomize with

  # use sigc as the template to randomize with
  sigc$chklen<-sigc$tend-sigc$tstart

  # randomization: for each chunk in sigc randomly choose an entry in pool_comb
  # that has length > the chklen
  # and tid does not exist in sigc 
  # randomly choose the start coord from 0 to v43_comb$tlen-sigc$chklen
  ck_chunks<-data.frame(tid=character(),
                        gid=character(),
                        name=character(),
                        tlen=numeric(),
                        tstart=numeric(),
                        tend=numeric())


  # remove the tids in sigc from the randomization pool
  pool_comb_fl<-pool_comb[(!pool_comb$tid %in% sigc$tid),]


  for (n in 1:nrow(sigc)) {
    
    temp.len<-sigc$chklen[n]
    
    pool_temp<-pool_comb_fl[which(pool_comb_fl$tlen>temp.len),]
    
    rand.row<-sample(c(1:nrow(pool_temp)),1)
    
    ck_row<-pool_temp[rand.row,]
    
    lendiff<-ck_row$tlen-temp.len
    
    
    if (ck_row$tid %in% ck_chunks$tid) {
      
      checkr<-ck_chunks[which(ck_chunks$tid==ck_row$tid),]
      
      excl<-checkr$tstart
      
      fran<-c(0:lendiff)
      
      valid_range <- setdiff(fran, excl)
      
      ck_row$tstart<-sample(valid_range,1)
      

    } else {
      
      ck_row$tstart<-sample(c(0:lendiff),1)
     
   }
    
    ck_row$tend<-ck_row$tstart+temp.len
    
    ck_chunks<-rbind(ck_chunks,ck_row)
    
  }


  # convert to bedfile
  ck_bedfile<-convert_transcript_to_genome(ck_chunks, exons)

  write.table(ck_bedfile,paste0('xist_sigchunk_',rp,'_ck_1.bed'),sep='\t', quote=F, col.names = F, row.names = F)


}




