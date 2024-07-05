# coparse chunk analysis
# set working directory to where the files are located

library(dplyr)
library(tidyr)

hrp4<-read.csv('v25_S4_v_XIST_repeats_pval.csv',header=T)
hrp5<-read.csv('v25_S5_v_XIST_repeats_pval.csv',header=T)

mrp4<-read.csv('vM10_S4_v_XIST_repeats_pval.csv',header=T)
mrp5<-read.csv('vM10_S5_v_XIST_repeats_pval.csv',header=T)

# convert the pval dataframe to sig chunk count dataframe for each repeats
# each row is a gene an columns are the sig chunk counts for each Xist repeats
labelchunk<-function(df) {
  
  chunknames<-strsplit(df$X,'|_', fixed=T)
  
  cl <- sapply(chunknames, length)
  print(sum(cl!=2))
  
  df$lncRNA<-sapply(chunknames,'[[',1)
  
  df$chunkstart<-0
  
  df$chunkstart[grepl('|_0',df$X,fixed=T)]<-1
  
  df$genenum<-cumsum(df$chunkstart)
  
  df$uniID<-paste0(df$lncRNA,'|_',df$genenum)
  
  count_sig<- function(x) sum(x < 0.05)
  
  siglnc <- aggregate(cbind(rA, rF, rB1, rB2, rD, rE) ~ uniID, data = df, FUN = count_sig)
  
  siglnc$totalcount<-rowSums(siglnc[,c(2:7)])
  
  sigtemp<-strsplit(siglnc$uniID,'|_',fixed=T)
  
  siglnc$lncRNA<-sapply(sigtemp,'[[',1)
  
  siglnc$lncRNA<-paste0(siglnc$lncRNA,'|')
  
  siglnc$order<-as.numeric(sapply(sigtemp,'[[',2))
  
  siglnc<-siglnc[order(siglnc$order,decreasing = FALSE),]
  
  return(siglnc)
}

hrp4_gene<-labelchunk(hrp4)
write.csv(hrp4_gene,'v25_S4_v_XIST_repeats_count.csv',row.names = F)

hrp5_gene<-labelchunk(hrp5)
write.csv(hrp5_gene,'v25_S5_v_XIST_repeats_count.csv',row.names = F)

mrp4_gene<-labelchunk(mrp4)
write.csv(mrp4_gene,'vM10_S4_v_XIST_repeats_count.csv',row.names = F)

mrp5_gene<-labelchunk(mrp5)
write.csv(mrp5_gene,'vM10_S5_v_XIST_repeats_count.csv',row.names = F)

# organize data to combine results together with previous S4 and S5 table

# load in human and mouse S4 compared with full XIST results
hxist_rval<-read.csv('v25_S4_vs_XIST_rval.csv',header=T)

hxist_pval<-read.csv('v25_S4_v_XIST_pval.csv',header=T)

hxist_adjpval<-read.csv('v25_S4_v_XIST_pval_bh.csv',header=T)


hxist_rval$X<-gsub('>','',hxist_rval$X)
sum(hxist_rval$X!=hxist_pval$X)
sum(hxist_rval$X!=hxist_adjpval$X)

colnames(hxist_rval)[2]<-'human_vs_XIST_rval'
colnames(hxist_pval)[2]<-'human_vs_XIST_pval'
colnames(hxist_adjpval)[2]<-'human_vs_XIST_adjpval'

colnames(hxist_pval)[1]<-'human_full_ID'


mxist_rval<-read.csv('vM10_S4_vs_XIST_rval.csv',header=T)

mxist_pval<-read.csv('vM10_S4_v_XIST_pval.csv',header=T)

mxist_adjpval<-read.csv('vM10_S4_v_XIST_pval_bh.csv',header=T)

mxist_rval$X<-gsub('>','',mxist_rval$X)
sum(mxist_rval$X!=mxist_pval$X)
sum(mxist_rval$X!=mxist_adjpval$X)


colnames(mxist_rval)[2]<-'mouse_vs_XIST_rval'
colnames(mxist_pval)[2]<-'mouse_vs_XIST_pval'
colnames(mxist_adjpval)[2]<-'mouse_vs_XIST_adjpval'

colnames(mxist_pval)[1]<-'mouse_full_ID'

# load in previous results and combine together

s4<-read.csv('coparse_tableS4_seekr.csv',header=T)

s4<-cbind(s4,hxist_pval$human_full_ID,hxist_rval$human_vs_XIST_rval,
          hxist_pval$human_vs_XIST_pval,hxist_adjpval$human_vs_XIST_adjpval,
          mxist_pval$mouse_full_ID,mxist_rval$mouse_vs_XIST_rval,
          mxist_pval$mouse_vs_XIST_pval,mxist_adjpval$mouse_vs_XIST_adjpval)


colnames(s4)<-c('human_transcript','human_tlen','mouse_transcript','mouse_tlen',
                'human_vs_mouse_rval','human_vs_mouse_pval','human_vs_mouse_adjpval',
                'human_full_ID','human_vs_XIST_rval','human_vs_XIST_pval','human_vs_XIST_adjpval',
                'mouse_full_ID','mouse_vs_XIST_rval','mouse_vs_XIST_pval','mouse_vs_XIST_adjpval')



hrp_count<-read.csv('v25_S4_v_XIST_repeats_count.csv',header=T)

sum(hrp_count$lncRNA!=hxist_rval$X)

hrp_count$uniID<-NULL
hrp_count$lncRNA<-NULL
hrp_count$order<-NULL

colnames(hrp_count)<-c('human_vs_rA_chunkcount','human_vs_rF_chunkcount','human_vs_rB1_chunkcount',
                       'human_vs_rB2_chunkcount','human_vs_rD_chunkcount','human_vs_rE_chunkcount',
                       'human_vs_XIST_total_chunkcount')

s4<-cbind(s4,hrp_count)


mrp_count<-read.csv('vM10_S4_v_XIST_repeats_count.csv',header=T)

sum(mrp_count$lncRNA!=mxist_rval$X)

mrp_count$uniID<-NULL
mrp_count$lncRNA<-NULL
mrp_count$order<-NULL


colnames(mrp_count)<-c('mouse_vs_rA_chunkcount','mouse_vs_rF_chunkcount','mouse_vs_rB1_chunkcount',
                       'mouse_vs_rB2_chunkcount','mouse_vs_rD_chunkcount','mouse_vs_rE_chunkcount',
                       'mouse_vs_XIST_total_chunkcount')


s4<-cbind(s4,mrp_count)

write.csv(s4,'coparse_tableS4_seekr_full.csv',row.names = F)


################################################################
# repeat the organization for S5 table

hxist_rval<-read.csv('v25_S5_vs_XIST_rval.csv',header=T)

hxist_pval<-read.csv('v25_S5_v_XIST_pval.csv',header=T)

hxist_adjpval<-read.csv('v25_S5_v_XIST_pval_bh.csv',header=T)


hxist_rval$X<-gsub('>','',hxist_rval$X)
sum(hxist_rval$X!=hxist_pval$X)
sum(hxist_rval$X!=hxist_adjpval$X)

colnames(hxist_rval)[2]<-'human_vs_XIST_rval'
colnames(hxist_pval)[2]<-'human_vs_XIST_pval'
colnames(hxist_adjpval)[2]<-'human_vs_XIST_adjpval'

colnames(hxist_pval)[1]<-'human_full_ID'


mxist_rval<-read.csv('vM10_S5_vs_XIST_rval.csv',header=T)

mxist_pval<-read.csv('vM10_S5_v_XIST_pval.csv',header=T)

mxist_adjpval<-read.csv('vM10_S5_v_XIST_pval_bh.csv',header=T)

mxist_rval$X<-gsub('>','',mxist_rval$X)
sum(mxist_rval$X!=mxist_pval$X)
sum(mxist_rval$X!=mxist_adjpval$X)


colnames(mxist_rval)[2]<-'mouse_vs_XIST_rval'
colnames(mxist_pval)[2]<-'mouse_vs_XIST_pval'
colnames(mxist_adjpval)[2]<-'mouse_vs_XIST_adjpval'

colnames(mxist_pval)[1]<-'mouse_full_ID'

s5<-read.csv('coparse_tableS5_seekr.csv',header=T)

s5<-cbind(s5,hxist_pval$human_full_ID,hxist_rval$human_vs_XIST_rval,
          hxist_pval$human_vs_XIST_pval,hxist_adjpval$human_vs_XIST_adjpval,
          mxist_pval$mouse_full_ID,mxist_rval$mouse_vs_XIST_rval,
          mxist_pval$mouse_vs_XIST_pval,mxist_adjpval$mouse_vs_XIST_adjpval)


colnames(s5)<-c('human_transcript','human_tlen','mouse_transcript','mouse_tlen',
                'human_vs_mouse_rval','human_vs_mouse_pval','human_vs_mouse_adjpval',
                'human_full_ID','human_vs_XIST_rval','human_vs_XIST_pval','human_vs_XIST_adjpval',
                'mouse_full_ID','mouse_vs_XIST_rval','mouse_vs_XIST_pval','mouse_vs_XIST_adjpval')



hrp_count<-read.csv('v25_S5_v_XIST_repeats_count.csv',header=T)

sum(hrp_count$lncRNA!=hxist_rval$X)

hrp_count$uniID<-NULL
hrp_count$lncRNA<-NULL
hrp_count$order<-NULL

colnames(hrp_count)<-c('human_vs_rA_chunkcount','human_vs_rF_chunkcount','human_vs_rB1_chunkcount',
                       'human_vs_rB2_chunkcount','human_vs_rD_chunkcount','human_vs_rE_chunkcount',
                       'human_vs_XIST_total_chunkcount')

s5<-cbind(s5,hrp_count)


mrp_count<-read.csv('vM10_S5_v_XIST_repeats_count.csv',header=T)

sum(mrp_count$lncRNA!=mxist_rval$X)

mrp_count$uniID<-NULL
mrp_count$lncRNA<-NULL
mrp_count$order<-NULL


colnames(mrp_count)<-c('mouse_vs_rA_chunkcount','mouse_vs_rF_chunkcount','mouse_vs_rB1_chunkcount',
                       'mouse_vs_rB2_chunkcount','mouse_vs_rD_chunkcount','mouse_vs_rE_chunkcount',
                       'mouse_vs_XIST_total_chunkcount')


s5<-cbind(s5,mrp_count)

write.csv(s5,'coparse_tableS5_seekr_full.csv',row.names = F)



