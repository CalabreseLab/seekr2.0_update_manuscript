# process multicov results 

# set workign directory to where the files are saved
library(dplyr)


# previously if a chunk overlaps with multiple exons
# it would be split into multiple entries according to the number of exons overlapped
# now we need to merge the counts for these split entries back to that chunk
# so that we can compare the sig chunk counts with its randomized control

merge_splits<-function(df) {
  
  df_merged <- df %>%
    group_by(name) %>%
    summarize(
      chrom = unique(chrom),
      chromStart = min(chromStart),
      chromEnd = max(chromEnd),
      score = unique(score),
      strand = unique(strand),
      RBM15_eclip = sum(RBM15_eclip),
      HNRNPK_eclip = sum(HNRNPK_eclip),
      MATR3_eclip = sum(MATR3_eclip),
      PTBP1_eclip = sum(PTBP1_eclip),
      HNRNPM_eclip = sum(HNRNPM_eclip),
      len = sum(len)
    ) %>%
    arrange(match(name, unique(df$name)))  # Preserving the original order
  df_merged<-as.data.frame(df_merged)
  return(df_merged)
}

# a separate function for merging splits in control as the colnames are different
merge_splits_ck<-function(df) {
  
  df_merged <- df %>%
    group_by(name_ck) %>%
    summarize(
      chrom_ck = unique(chrom_ck),
      chromStart_ck = min(chromStart_ck),
      chromEnd_ck = max(chromEnd_ck),
      score_ck = unique(score_ck),
      strand_ck = unique(strand_ck),
      RBM15_ck = sum(RBM15_ck),
      HNRNPK_ck = sum(HNRNPK_ck),
      MATR3_ck = sum(MATR3_ck),
      PTBP1_ck = sum(PTBP1_ck),
      HNRNPM_ck = sum(HNRNPM_ck),
      len_ck = sum(len_ck)
    ) %>%
    arrange(match(name_ck, unique(df$name_ck)))  # Preserving the original order
  df_merged<-as.data.frame(df_merged)
  return(df_merged)
}

# read in the multicov output for sig chunks and rename the columns
sigc_a<-read.table('chunk_rA_multicov.out',sep='\t')
colnames(sigc_a)<-c('chrom','chromStart','chromEnd','name','score','strand','RBM15_eclip','HNRNPK_eclip','MATR3_eclip','PTBP1_eclip','HNRNPM_eclip')
# add a length column
sigc_a$len<-sigc_a$chromEnd-sigc_a$chromStart+1
sigc_a<-merge_splits(sigc_a)

# read in the control file
ck_a<-read.table('control_rA_multicov.out',sep='\t')
colnames(ck_a)<-c('chrom_ck','chromStart_ck','chromEnd_ck','name_ck','score_ck','strand_ck','RBM15_ck','HNRNPK_ck','MATR3_ck','PTBP1_ck','HNRNPM_ck')
ck_a$len_ck<-ck_a$chromEnd_ck-ck_a$chromStart_ck+1
ck_a<-merge_splits_ck(ck_a)

# check if the length of sigc is the same as control -- should be all the same 
# expect outcome: 0
sum(sigc_a$len!=ck_a$len_ck)

# combine sigc and ck together and save as a csv file
rA<-cbind(sigc_a,ck_a)

write.csv(rA,'chunk_multicov_rA_comb.csv',row.names = F)



sigc_f<-read.table('chunk_rF_multicov.out',sep='\t')
colnames(sigc_f)<-c('chrom','chromStart','chromEnd','name','score','strand','RBM15_eclip','HNRNPK_eclip','MATR3_eclip','PTBP1_eclip','HNRNPM_eclip')
sigc_f$len<-sigc_f$chromEnd-sigc_f$chromStart+1
sigc_f<-merge_splits(sigc_f)


ck_f<-read.table('control_rF_multicov.out',sep='\t')
colnames(ck_f)<-c('chrom_ck','chromStart_ck','chromEnd_ck','name_ck','score_ck','strand_ck','RBM15_ck','HNRNPK_ck','MATR3_ck','PTBP1_ck','HNRNPM_ck')
ck_f$len_ck<-ck_f$chromEnd_ck-ck_f$chromStart_ck+1
ck_f<-merge_splits_ck(ck_f)

sum(sigc_f$len!=ck_f$len_ck)

rF<-cbind(sigc_f,ck_f)

write.csv(rF,'chunk_multicov_rF_comb.csv',row.names = F)



sigc_b1<-read.table('chunk_rB1_multicov.out',sep='\t')
colnames(sigc_b1)<-c('chrom','chromStart','chromEnd','name','score','strand','RBM15_eclip','HNRNPK_eclip','MATR3_eclip','PTBP1_eclip','HNRNPM_eclip')
sigc_b1$len<-sigc_b1$chromEnd-sigc_b1$chromStart+1
sigc_b1<-merge_splits(sigc_b1)


ck_b1<-read.table('control_rB1_multicov.out',sep='\t')
colnames(ck_b1)<-c('chrom_ck','chromStart_ck','chromEnd_ck','name_ck','score_ck','strand_ck','RBM15_ck','HNRNPK_ck','MATR3_ck','PTBP1_ck','HNRNPM_ck')
ck_b1$len_ck<-ck_b1$chromEnd_ck-ck_b1$chromStart_ck+1
ck_b1<-merge_splits_ck(ck_b1)

sum(sigc_b1$len!=ck_b1$len_ck)


rB1<-cbind(sigc_b1,ck_b1)

write.csv(rB1,'chunk_multicov_rB1_comb.csv',row.names = F)


sigc_b2<-read.table('chunk_rB2_multicov.out',sep='\t')
colnames(sigc_b2)<-c('chrom','chromStart','chromEnd','name','score','strand','RBM15_eclip','HNRNPK_eclip','MATR3_eclip','PTBP1_eclip','HNRNPM_eclip')
sigc_b2$len<-sigc_b2$chromEnd-sigc_b2$chromStart+1
sigc_b2<-merge_splits(sigc_b2)

ck_b2<-read.table('control_rB2_multicov.out',sep='\t')
colnames(ck_b2)<-c('chrom_ck','chromStart_ck','chromEnd_ck','name_ck','score_ck','strand_ck','RBM15_ck','HNRNPK_ck','MATR3_ck','PTBP1_ck','HNRNPM_ck')
ck_b2$len_ck<-ck_b2$chromEnd_ck-ck_b2$chromStart_ck+1
ck_b2<-merge_splits_ck(ck_b2)

sum(sigc_b2$len!=ck_b2$len_ck)

rB2<-cbind(sigc_b2,ck_b2)

write.csv(rB2,'chunk_multicov_rB2_comb.csv',row.names = F)


sigc_d<-read.table('chunk_rD_multicov.out',sep='\t')
colnames(sigc_d)<-c('chrom','chromStart','chromEnd','name','score','strand','RBM15_eclip','HNRNPK_eclip','MATR3_eclip','PTBP1_eclip','HNRNPM_eclip')
sigc_d$len<-sigc_d$chromEnd-sigc_d$chromStart+1
sigc_d<-merge_splits(sigc_d)

ck_d<-read.table('control_rD_multicov.out',sep='\t')
colnames(ck_d)<-c('chrom_ck','chromStart_ck','chromEnd_ck','name_ck','score_ck','strand_ck','RBM15_ck','HNRNPK_ck','MATR3_ck','PTBP1_ck','HNRNPM_ck')
ck_d$len_ck<-ck_d$chromEnd_ck-ck_d$chromStart_ck+1
ck_d<-merge_splits_ck(ck_d)

sum(sigc_d$len!=ck_d$len_ck)


rD<-cbind(sigc_d,ck_d)

write.csv(rD,'chunk_multicov_rD_comb.csv',row.names = F)


sigc_e<-read.table('chunk_rE_multicov.out',sep='\t')
colnames(sigc_e)<-c('chrom','chromStart','chromEnd','name','score','strand','RBM15_eclip','HNRNPK_eclip','MATR3_eclip','PTBP1_eclip','HNRNPM_eclip')
sigc_e$len<-sigc_e$chromEnd-sigc_e$chromStart+1
sigc_e<-merge_splits(sigc_e)

ck_e<-read.table('control_rE_multicov.out',sep='\t')
colnames(ck_e)<-c('chrom_ck','chromStart_ck','chromEnd_ck','name_ck','score_ck','strand_ck','RBM15_ck','HNRNPK_ck','MATR3_ck','PTBP1_ck','HNRNPM_ck')
ck_e$len_ck<-ck_e$chromEnd_ck-ck_e$chromStart_ck+1
ck_e<-merge_splits_ck(ck_e)

sum(sigc_e$len!=ck_e$len_ck)


rE<-cbind(sigc_e,ck_e)

write.csv(rE,'chunk_multicov_rE_comb.csv',row.names = F)


###################################################################
# plotting


library(ggplot2)
library(tidyr)

# can be used to potentially filter for counts
custom.filter<-function(x,y){
  
  temp<-as.data.frame(cbind(x,y))
  #temp<-temp[which(temp[,1]>0 | temp[,2]>0),]  
  return(temp)
}


sigcstat<-as.data.frame(matrix(nrow=6,ncol=5))
rownames(sigcstat)<-c('rA','rB1','rB2','rD','rE','rF')
colnames(sigcstat)<-c('RBM15','HNRNPK','MATR3','PTBP1','HNRNPM')

sigccount<-as.data.frame(matrix(nrow=6,ncol=5))
rownames(sigccount)<-c('rA','rB1','rB2','rD','rE','rF')
colnames(sigccount)<-c('RBM15','HNRNPK','MATR3','PTBP1','HNRNPM')

eclip.vec<-c('RBM15_eclip','HNRNPK_eclip','MATR3_eclip','PTBP1_eclip','HNRNPM_eclip')
ck.vec<-c('RBM15_ck','HNRNPK_ck','MATR3_ck','PTBP1_ck','HNRNPM_ck')
repeat.vec<-c('rA','rB1','rB2','rD','rE','rF')

# remove variable if already exists
if exists(comb_long) {
  rm(comb_long)
}
if exists(comb_diff) {
  rm(comb_diff)
}


# organize stats for plotting
# loop through all repeats and RBPs
for (rp in 1:length(repeat.vec)) {
  # Use get() to assign the corresponding dataframe to df
  rpname<-repeat.vec[rp]
  df <- get(rpname)
  
  for (n in 1:length(eclip.vec)) {
    
    eclip<-eclip.vec[n]
    ck<-ck.vec[n]
    
    temp<-custom.filter(df[eclip],df[ck])
    # wilcox paired test to compare counts of sigc vs ck
    sigcstat[rp,n]<-wilcox.test(temp[,1],temp[,2],alternative = 'greater',paired=T)$p.value
    sigccount[rp,n]<-nrow(temp)
    temp$repts<-rpname
    temp$bedid<-row.names(temp)
    # convert to long format for plotting purpose
    temp_long <- pivot_longer(temp, cols = -c(repts,bedid),names_to = "exp", values_to = "count")
    temp_long<-as.data.frame(temp_long)
    if (exists('comb_long')) {
      comb_long<-rbind(comb_long,temp_long)
    } else {
      comb_long<-temp_long
    }
    
    
  }
  
}

# adjust p value
sigcstatadj<-p.adjust(as.vector(t(sigcstat)) ,method='bonferroni')
sigcstatadj<-matrix(sigcstatadj, nrow = 6, ncol = 5, byrow = TRUE)  
rownames(sigcstatadj)<-c('rA','rB1','rB2','rD','rE','rF')
colnames(sigcstatadj)<-c('RBM15','HNRNPK','MATR3','PTBP1','HNRNPM')



write.csv(sigcstatadj,'chunk_multicov_comb_paired_adjpval.csv')
write.csv(sigccount,'chunk_multicov_comb_paired_ncount.csv')


comb_long$exp<-factor(comb_long$exp,levels=c('RBM15_eclip','RBM15_ck','HNRNPK_eclip','HNRNPK_ck',
                                             'MATR3_eclip','MATR3_ck','PTBP1_eclip','PTBP1_ck',
                                             'HNRNPM_eclip','HNRNPM_ck'))

# plot y axis on the log scale -- log10(y+1)
sigcc<-ggplot(comb_long, aes(x = exp, y = (count+1), fill=exp)) +
  geom_boxplot(color='black',aes(fill=exp)) +
  scale_y_log10()+
  theme(
    panel.background=element_rect(fill='white'),
    panel.grid.major=element_line(color='grey',size=0.3),
    legend.position = 'none',
    legend.key.height = unit(1,'line'),
    legend.key.width  = unit(1,'line'),
    legend.title = element_text(size=32),
    legend.text=element_text(size=26),
    legend.text.align=0,
    axis.text.x=element_text(size=26,color='black',angle=45,hjust=1),
    axis.text.y=element_text(size=26,color='black'),
    axis.title.x=element_text(size=32),
    axis.title.y=element_text(size=32,angle=90),
    strip.text = element_text(size = 32)
  ) +
  labs(x = "eclip experiments", y = "Count")+
  scale_fill_brewer(palette = 'Paired')+
  facet_wrap(~repts, ncol = 3, scales = "free")

pdf('chunk_boxplot_all_count_2row.pdf',width=32,height=18)
plot(sigcc)
dev.off()
