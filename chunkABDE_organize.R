# organize the spreadsheet to see how each lncRNA is similar to ABDEF repeats
# setwd: set the working directory to the folder that contains the file v43_chunks_v_ABDEF_pvals.csv

chunk<-read.csv('v43_chunks_v_ABDEF_pvals.csv', header=T, row.names = 1)

chunknames<-row.names(chunk)

chunknames<-strsplit(chunknames,'|_', fixed=T)

cl <- sapply(chunknames, length)
sum(cl==2)

chunk$lncRNA<-sapply(chunknames,'[[',1)

count_sig<- function(x) sum(x < 0.05)

siglnc <- aggregate(cbind(rA, rF, rB1, rB2, rD, rE) ~ lncRNA, data = chunk, FUN = count_sig)

siglnc$countsums<-rowSums(siglnc[,c('rA','rF','rB1','rB2','rD','rE')])

siglnc<-siglnc[order(siglnc$countsums,decreasing = T),]

write.csv(siglnc,'v43_lncRNA_ABDEF_count.csv',row.names = F)

