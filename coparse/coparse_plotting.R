# plots for coparse data
# set working directory to the folder where all the files are saved


library(ggplot2)
library(Biostrings) # needed to be loaded first, before rtracklayer
library(rtracklayer)
library(ggrepel)

# read in data
hh <- read.csv('human_vs_human_6mer_rvals_1m.csv',header=T)

colnames(hh)<-c('rowname','rval')

s4r <- read.csv('s4_diag_rvals.csv',header=T)
s5r <- read.csv('s5_diag_rvals.csv',header=T)

s4p <- read.csv('s4_diag_pvals.csv',header=T)
s5p <- read.csv('s5_diag_pvals.csv',header=T)

s4adjp <- read.csv('s4_diag_bh_adjpvals.csv',header=T)
s5adjp <- read.csv('s5_diag_bh_adjpvals.csv',header=T)

s4df<-read.csv('coparse_tableS4.csv',header=T)
s5df<-read.csv('coparse_tableS5a.csv',header=T)

# only keep the human/mouse transcript columns for simplicity
s4df<-s4df[,c('human.transcript','mouse.transcript')]
s5df<-s5df[,c('human.transcript','mouse.transcript')]

# combine the transcript and r p adjp values all together into one dataframe
s4df<-cbind(s4df,s4r,s4p,s4adjp)
colnames(s4df)<-c('human_transcript','mouse_transcript','rval','pval','adjpval')


s5df<-cbind(s5df,s5r,s5p,s5adjp)
colnames(s5df)<-c('human_transcript','mouse_transcript','rval','pval','adjpval')


#########################################################
# wilcox test

wilcox.test(hh$rval,s4df$rval,alternative = 'less')

wilcox.test(hh$rval,s5df$rval,alternative = 'less')

##########################################################
# add transcript length to the table
s4h_seq<-readDNAStringSet('v25_S4.fa')
s4m_seq<-readDNAStringSet('vM10_S4.fa')

s5h_seq<-readDNAStringSet('v25_S5.fa')
s5m_seq<-readDNAStringSet('vM10_S5.fa')

s4h_name<-names(s4h_seq)
s4m_name<-names(s4m_seq)

s5h_name<-names(s5h_seq)
s5m_name<-names(s5m_seq)

rm(s4h_seq)
rm(s4m_seq)
rm(s5h_seq)
rm(s5m_seq)

# organize the data and convert the fasta header into transcriptID, transcript name and transcript length
org_data<-function(vec) {
  
  vec<-strsplit(vec,'|',fixed=T)
  print(sum(lengths(vec)!=7))
  
  transcript_ID<-sapply(vec,'[[',1)
  transcript_ID<-strsplit(transcript_ID,'.',fixed=T)
  transcript_ID<-sapply(transcript_ID,'[[',1)
  
  tname<-sapply(vec,'[[',5)
  tlen<-sapply(vec,'[[',7)
  
  df<-as.data.frame(cbind(transcript_ID,tname,tlen))
  df$tlen<-as.numeric(df$tlen)
  
  return(df)
  
}

s4h_data<-org_data(s4h_name)
s4m_data<-org_data(s4m_name)

s5h_data<-org_data(s5h_name)
s5m_data<-org_data(s5m_name)

s4df$human_tlen<-NA
s4df$mouse_tlen<-NA

# adding in the transcript length for human and mouse for S4
for (n in 1:nrow(s4df)) {
  htemp<-s4df$human_transcript[n]
  mtemp<-s4df$mouse_transcript[n]
  
  s4df$human_tlen[n]<-unique(s4h_data$tlen[which(s4h_data$transcript_ID==htemp)])
  
  s4df$mouse_tlen[n]<-unique(s4m_data$tlen[which(s4m_data$transcript_ID==mtemp)])
  
}


s5df$human_tlen<-NA
s5df$mouse_tlen<-NA

# adding in the transcript length for human and mouse for S5
for (n in 1:nrow(s5df)) {
  htemp<-s5df$human_transcript[n]
  mtemp<-s5df$mouse_transcript[n]
  
  s5df$human_tlen[n]<-unique(s5h_data$tlen[which(s5h_data$transcript_ID==htemp)])
  
  s5df$mouse_tlen[n]<-unique(s5m_data$tlen[which(s5m_data$transcript_ID==mtemp)])
  
}

# re-organize the dataframe to have the orders show in colnames below
s4df<-s4df[,c(1,6,2,7,3,4,5)]
# colnames(s4df)<-c('human_transcript','human_tlen','mouse_transcript','mouse_tlen','rval','pval','adjpval')
s5df<-s5df[,c(1,6,2,7,3,4,5)]
# colnames(s5df)<-c('human_transcript','human_tlen','mouse_transcript','mouse_tlen','rval','pval','adjpval')


write.csv(s4df,'coparse_tableS4_seekr.csv',row.names = F)
write.csv(s5df,'coparse_tableS5_seekr.csv',row.names = F)

########################################################
# read in the highly conserved lncRNAs
hc<-read.table('highly_conserved_lncRNAS.txt',sep=' ')
hc$V1<-NULL
colnames(hc)<-c('human_name','mouse_name')

# read in the fasta files to extract sequences
v25lnc<-readDNAStringSet('/Users/shuang/Downloads/gencode.v25.lncRNA_transcripts.fa')
vM10lnc<-readDNAStringSet('/Users/shuang/Downloads/gencode.vM10.lncRNA_transcripts.fa')

# parse the fasta header to match human/mouse_name
v25lnc_name<-names(v25lnc)
vM10lnc_name<-names(vM10lnc)

v25lnc_name<-strsplit(v25lnc_name,'|',fixed=T)
vM10lnc_name<-strsplit(vM10lnc_name,'|',fixed=T)

v25lnc_name<-sapply(v25lnc_name,'[[',5)
vM10lnc_name<-sapply(vM10lnc_name,'[[',5)

hc$hidx<-NA
hc$midx<-NA

hc<-hc[c(1:13,15,14),]

for (n in 1:(nrow(hc)-1)) {
  
  hc$hidx[n]<-which(v25lnc_name==hc$human_name[n])
  hc$midx[n]<-which(vM10lnc_name==hc$mouse_name[n])
  
}

# save the sequences
human_hc <- v25lnc[hc$hidx[1:14]]
writeXStringSet(human_hc, filepath = "human_highly_conserved.fa")

mouse_hc<-vM10lnc[hc$midx[1:14]]
writeXStringSet(mouse_hc, filepath = "mouse_highly_conserved.fa")

# get TUNAR sequence from the whole transcript fasta file and append to the end

hcrvals<-read.csv('highly_conserved_rvals.csv',header=T)

hc$rval<-hcrvals$X0

hc$label<-gsub('-001','',hc$human_name)

hcall<-hc

ss<-c('HOTAIR','HOTTIP','TUG1','H19','NEAT1','XIST','MALAT1','KCNQ1OT1')

hc<-hc[(hc$label %in% ss),]

########################################################
# set coeff as the ratio between the first y axis and the second y axis
# this can be calculate maximum counts for histogram divided by the max density
# coeff<-(max_counts / max_density)

coeff<-150

# set histogram color based on adjpval
s4df$fillcol<-'white'
s4df$fillcol[which(s4df$adjpval<0.05)]<-'black'
# find the r val on the cutoff
min(s4df$rval[which(s4df$fillcol=='black')]) # 0.1575171
max(s4df$rval[which(s4df$fillcol=='white')]) # 0.1574389
#breaks has to be on 0.1575

hbreaks<-seq(from=-1.0425, to=1.0575, by=0.02)


# Create the plot
p<-ggplot() +
  geom_density(data = hh, aes(x = rval, y=after_stat(density)), fill="darkgrey",color='black', alpha=0.8) +
  geom_histogram(data = s4df, aes(x = rval, y = after_stat(count)/coeff, fill=fillcol,color=fillcol), breaks=hbreaks, alpha=0.5) +
  scale_y_continuous(name = "Density",breaks=seq(from=0, to=12, by=2), sec.axis = sec_axis(~ . * coeff,name = "Counts",breaks = seq(0, 2400, by = 300))) +
  labs(x = "SEEKR R Values",y = "Density") +
  coord_cartesian(xlim=c(-0.2,0.7))+
  scale_x_continuous(breaks=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7))+
  #scale_y_continuous(breaks=seq(from=0, to=12, by=2))+
  scale_color_manual(values=c('#d7191c','#2c7bb6'))+
  scale_fill_manual(values=c('#fdae61','#abd9e9'))+
  theme(plot.title=element_text(size=22),
        panel.background=element_rect(fill='white'),
        #plot.margin = margin(2, 2, 2, 2, "pt"),
        # panel.grid.major=element_line(color='grey',linewidth =0.3),
        panel.grid.major=element_blank(),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title.y.right = element_text(size=20),
        axis.text.y.right = element_text(size=20))

p<-p + geom_vline(data = hc, aes(xintercept = rval), color = "#1b9e77", linetype = "solid",linewidth=0.5)+
  geom_text_repel(data = hc, aes(x = rval, y = 11, label = label), 
                  size = 4, nudge_x = 0.01, direction = "y", hjust = -0.1, vjust = -0.5,
                  segment.color = "darkgray",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid")

  

ggsave("S4_dens_hist_adjpval.pdf", plot = p, width = 6.5, height = 3, units = "in")


##################################
# set histogram color based on pval
coeff<-150

s4df$fillcol<-'white'
s4df$fillcol[which(s4df$pval<0.05)]<-'black'
# find the r val on the cutoff
min(s4df$rval[which(s4df$fillcol=='black')]) # 0.07072676
max(s4df$rval[which(s4df$fillcol=='white')]) # 0.07058748
#breaks has to be on 0.0706

hbreaks<-seq(from=-1.0294, to=1.0706, by=0.02)


# Create the plot
p<-ggplot() +
  geom_density(data = hh, aes(x = rval, y=after_stat(density)), fill="darkgrey",color='black', alpha=0.8) +
  geom_histogram(data = s4df, aes(x = rval, y = after_stat(count)/coeff, fill=fillcol,color=fillcol), breaks=hbreaks, alpha=0.5) +
  scale_y_continuous(name = "Density",breaks=seq(from=0, to=12, by=2), sec.axis = sec_axis(~ . * coeff,name = "Counts",breaks = seq(0, 2400, by = 300))) +
  labs(x = "SEEKR R Values",y = "Density") +
  coord_cartesian(xlim=c(-0.2,0.7))+
  scale_x_continuous(breaks=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7))+
  #scale_y_continuous(breaks=seq(from=0, to=12, by=2))+
  scale_color_manual(values=c('#d7191c','#2c7bb6'))+
  scale_fill_manual(values=c('#fdae61','#abd9e9'))+
  theme(plot.title=element_text(size=22),
        panel.background=element_rect(fill='white'),
        #plot.margin = margin(2, 2, 2, 2, "pt"),
        #panel.grid.major=element_line(color='grey',linewidth =0.3),
        panel.grid.major=element_blank(),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title.y.right = element_text(size=20),
        axis.text.y.right = element_text(size=20))

p<-p + geom_vline(data = hc, aes(xintercept = rval), color = "#1b9e77", linetype = "solid",linewidth=0.5)+
  geom_text_repel(data = hc, aes(x = rval, y = 11, label = label), 
                  size = 4, nudge_x = 0.01, direction = "y", hjust = -0.1, vjust = -0.5,
                  segment.color = "darkgray",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid")

ggsave("S4_dens_hist_pval.pdf", plot = p, width = 6.5, height = 3, units = "in")


#############################
# S5
#############################


# set coeff as the ratio between the first y axis and the second y axis
# this can be calculate maximum counts for histogram divided by the max density
# coeff<-(max_counts / max_density)
coeff<-10

# set histogram color based on adjpval
s5df$fillcol<-'white'
s5df$fillcol[which(s5df$adjpval<0.05)]<-'black'

# find the r val on the cutoff
min(s5df$rval[which(s5df$fillcol=='black')]) # 0.1413217
max(s5df$rval[which(s5df$fillcol=='white')]) # 0.1345429

# breaks has to be on 0.14
hbreaks<-seq(from=-1.06, to=1.04, by=0.02)


# Create the plot
p<-ggplot() +
  geom_density(data = hh, aes(x = rval, y=after_stat(density)), fill="darkgrey",color='black', alpha=0.8) +
  geom_histogram(data = s5df, aes(x = rval, y = after_stat(count)/coeff, fill=fillcol,color=fillcol), breaks=hbreaks, alpha=0.5) +
  scale_y_continuous(name = "Density",breaks=seq(from=0, to=12, by=2), sec.axis = sec_axis(~ . * coeff,name = "Counts",breaks = seq(0, 120, by = 20))) +
  labs(x = "SEEKR R Values",y = "Density") +
  coord_cartesian(xlim=c(-0.2,0.7))+
  scale_x_continuous(breaks=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7))+
  #scale_y_continuous(breaks=seq(from=0, to=12, by=2))+
  scale_color_manual(values=c('#d7191c','#2c7bb6'))+
  scale_fill_manual(values=c('#fdae61','#abd9e9'))+
  theme(plot.title=element_text(size=22),
        panel.background=element_rect(fill='white'),
        #plot.margin = margin(2, 2, 2, 2, "pt"),
        panel.grid.major=element_blank(),
        #panel.grid.major=element_line(color='grey',linewidth =0.3),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title.y.right = element_text(size=20),
        axis.text.y.right = element_text(size=20))

p<-p + geom_vline(data = hc, aes(xintercept = rval), color = "#1b9e77", linetype = "solid",linewidth=0.5)+
  geom_text_repel(data = hc, aes(x = rval, y = 11, label = label), 
                  size = 4, nudge_x = 0.01, direction = "y", hjust = -0.1, vjust = -0.5,
                  segment.color = "darkgray",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid")

ggsave("S5_dens_hist_adjpval.pdf", plot = p, width = 6.5, height = 3, units = "in")


#######################################
coeff<-10

# set histogram color based on adjpval
s5df$fillcol<-'white'
s5df$fillcol[which(s5df$pval<0.05)]<-'black'


# find the r val on the cutoff
min(s5df$rval[which(s5df$fillcol=='black')]) # 0.07088741
max(s5df$rval[which(s5df$fillcol=='white')]) # 0.0695067

# breaks has to be on 0.07
hbreaks<-seq(from=-1.03, to=1.07, by=0.02)


# Create the plot
p<-ggplot() +
  geom_density(data = hh, aes(x = rval, y=after_stat(density)), fill="darkgrey",color='black', alpha=0.8) +
  geom_histogram(data = s5df, aes(x = rval, y = after_stat(count)/coeff, fill=fillcol,color=fillcol), breaks=hbreaks, alpha=0.5) +
  scale_y_continuous(name = "Density",breaks=seq(from=0, to=12, by=2), sec.axis = sec_axis(~ . * coeff,name = "Counts",breaks = seq(0, 120, by = 20))) +
  labs(x = "SEEKR R Values",y = "Density") +
  coord_cartesian(xlim=c(-0.2,0.7))+
  scale_x_continuous(breaks=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7))+
  #scale_y_continuous(breaks=seq(from=0, to=12, by=2))+
  scale_color_manual(values=c('#d7191c','#2c7bb6'))+
  scale_fill_manual(values=c('#fdae61','#abd9e9'))+
  theme(plot.title=element_text(size=22),
        panel.background=element_rect(fill='white'),
        #plot.margin = margin(2, 2, 2, 2, "pt"),
        # panel.grid.major=element_line(color='grey',linewidth =0.3),
        panel.grid.major=element_blank(),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title.y.right = element_text(size=20),
        axis.text.y.right = element_text(size=20))


p<-p + geom_vline(data = hc, aes(xintercept = rval), color = "#1b9e77", linetype = "solid",linewidth=0.5)+
  geom_text_repel(data = hc, aes(x = rval, y = 11, label = label), 
                  size = 4, nudge_x = 0.01, direction = "y", hjust = -0.1, vjust = -0.5,
                  segment.color = "darkgray",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid")

ggsave("S5_dens_hist_pval.pdf", plot = p, width = 6.5, height = 3, units = "in")

