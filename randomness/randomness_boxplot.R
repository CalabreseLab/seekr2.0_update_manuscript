########## plot the all vs common10

library(ggplot2)

plotdf<-read.csv('allvscommon10.csv',header=T)

plotdf$test<-factor(plotdf$test,levels=c('common10','all'))

plotdf$model<-factor(plotdf$model,levels=c('lognorm','johnsonsu','genhyperbolic'))

p<-ggplot(data=plotdf,aes(x=test,y=r95)) +
  labs(x = "Test",y = "r value @ p=0.05") +
  scale_y_continuous(limits = c(0.07,0.08), breaks=seq(from=0.07, to=0.08, by=0.005)) +
  theme(plot.title=element_blank(),
        panel.background=element_rect(fill='white'),
        panel.grid.major=element_line(color='grey',linewidth =0.3),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.title=element_text(size=16),
        legend.position='top',
        legend.key.height = unit(1,'line'),
        legend.key.width  = unit(0.5,'line'),
        legend.key=element_rect(fill='transparent'),
        legend.text=element_text(size=16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16, angle = 45, hjust = 0.5),
        axis.title.y.right = element_text(size=16),
        axis.text.y.right = element_text(size=16)) +
  geom_boxplot(fill='grey',color='#708090',alpha=0.3,outlier.shape = NA,width = 0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "#708090") +
  geom_point(aes(x=test,y=r95,color=test,shape=model), 
             position = position_jitter(width=0.25,height=0), size = 3,alpha=0.7) +
  scale_color_brewer(palette = 'Set1',guide = 'none')+
  scale_shape_manual(values = c(16,17,15))+
  guides(shape=guide_legend(nrow=3))

ggsave("allvs10_r95.pdf", plot = p, width = 3, height = 3, units = "in")



p<-ggplot(data=plotdf,aes(x=test,y=psigcount)) +
  labs(x = "Test",y = "Sig Gene Count") +
  scale_y_continuous(limits = c(2700,3100), breaks=seq(from=2700, to=3100, by=200)) +
  theme(plot.title=element_blank(),
        panel.background=element_rect(fill='white'),
        panel.grid.major=element_line(color='grey',linewidth =0.3),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.title=element_text(size=16),
        legend.position='top',
        legend.key.height = unit(1,'line'),
        legend.key.width  = unit(0.5,'line'),
        legend.key=element_rect(fill='transparent'),
        legend.text=element_text(size=16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16, angle = 45, hjust = 0.5),
        axis.title.y.right = element_text(size=16),
        axis.text.y.right = element_text(size=16)) +
  geom_boxplot(fill='grey',color='#708090',alpha=0.3,outlier.shape = NA,width = 0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "#708090") +
  geom_point(aes(x=test,y=psigcount,color=test,shape=model), 
             position = position_jitter(width=0.25,height=0), size = 3,alpha=0.7) +
  scale_color_brewer(palette = 'Set1',guide = 'none')+
  scale_shape_manual(values = c(16,17,15))+
  guides(shape=guide_legend(nrow=3))

ggsave("allvs10_pcount.pdf", plot = p, width = 3, height = 3, units = "in")



# plot the randomness box plot
# set working directory to where the files are located

library(ggplot2)

# load in the p val and adjusted p val file names
files<-dir('./')

file_10k_p<-files[grepl("^v43vXIST_10k_(10|[1-9])\\.csv$",files)]

file_100k_p<-files[grepl("^v43vXIST_100k_(10|[1-9])\\.csv$",files)]

file_1M_p<-files[grepl("^v43vXIST_1M_(10|[1-9])\\.csv$",files)]

file_full_p<-files[grepl("^v43vXIST_full_(10|[1-9])\\.csv$",files)]

file_10k_adjp<-files[grepl("^v43vXIST_10k_bh_(10|[1-9])\\.csv$",files)]

file_100k_adjp<-files[grepl("^v43vXIST_100k_bh_(10|[1-9])\\.csv$",files)]

file_1M_adjp<-files[grepl("^v43vXIST_1M_bh_(10|[1-9])\\.csv$",files)]

file_full_adjp<-files[grepl("^v43vXIST_full_bh_(10|[1-9])\\.csv$",files)]


# count the significant gene number
count_sig<-function(filename) {
  
  df<-read.csv(filename,header=T)
  
  count<-sum(df$ENST00000429829.6_XIST<0.05)
  
  return(count)
}


# organize data for plotting
subsetsize<-rep(c('10k','100k','1M','full'),each=10)
plotdf<-as.data.frame(subsetsize)
plotdf$count<-NA

filecomb_p<-c(file_10k_p,file_100k_p,file_1M_p,file_full_p)

# add the sig count based on p val
for (n in 1:length(filecomb_p)) {
  
  f<-filecomb_p[n]
  
  plotdf$count[n]<-count_sig(f)
}

# set the levels for the subset size 
plotdf$subsetsize<-factor(plotdf$subsetsize,levels=c('10k','100k','1M','full'))

pvalcomb<-plotdf

# get the data for adjusted p val
subsetsize<-rep(c('10k','100k','1M','full'),each=10)
plotdf<-as.data.frame(subsetsize)
plotdf$count<-NA

filecomb_adjp<-c(file_10k_adjp,file_100k_adjp,file_1M_adjp,file_full_adjp)

for (n in 1:length(filecomb_adjp)) {
  
  f<-filecomb_adjp[n]
  
  plotdf$count[n]<-count_sig(f)
}

plotdf$subsetsize<-factor(plotdf$subsetsize,levels=c('10k','100k','1M','full'))

sum(plotdf$subsetsize!=pvalcomb$subsetsize)

colnames(pvalcomb)<-c('subsetsize','pvalcount')

# combine p val and adjusted p val together
pvalcomb<-cbind(pvalcomb,plotdf$count)
colnames(pvalcomb)[3]<-'adjpvalcount'

# add in the best model (first one) used for each simulation 

# load in the mode fitting file names
files<-dir('./')

model_10k<-files[grepl("^v43_10k_(10|[1-9])\\.csv$",files)]
model_100k<-files[grepl("^v43_100k_(10|[1-9])\\.csv$",files)]
model_1M<-files[grepl("^v43_1M_(10|[1-9])\\.csv$",files)]
model_full<-files[grepl("^v43_full_(10|[1-9])\\.csv$",files)]

# get the best fit model name
get_bf<-function(filename) {
  
  df<-read.csv(filename,header=T)
  
  bf<-df$distribution_name[1]
  
  return(bf)
}

pvalcomb$model<-NA

filecomb_m<-c(model_10k,model_100k,model_1M,model_full)

for (n in 1:length(filecomb_m)) {
  
  f<-filecomb_m[n]
  
  pvalcomb$model[n]<-get_bf(f)
}

write.csv(pvalcomb,'randomnesscomb.csv',row.names = F)

pvalcomb$subsetsize<-factor(pvalcomb$subsetsize,levels=c('10k','100k','1M','full'))

pvalcomb$model<-factor(pvalcomb$model,levels=c('lognorm','chi2','rayleigh'))

p<-ggplot(data=pvalcomb,aes(x=subsetsize,y=pvalcount)) +
  labs(x = "Subset Size",y = "Sig Gene Count") +
  scale_y_continuous(limits = c(2500,3000), breaks=seq(from=2500, to=3000, by=250)) +
  theme(plot.title=element_blank(),
        panel.background=element_rect(fill='white'),
        panel.grid.major=element_line(color='grey',linewidth =0.3),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.title=element_text(size=10),
        legend.position='top',
        legend.key.height = unit(1,'line'),
        legend.key.width  = unit(0.5,'line'),
        legend.key=element_rect(fill='transparent'),
        legend.text=element_text(size=10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.y.right = element_text(size=10),
        axis.text.y.right = element_text(size=10)) +
  geom_boxplot(fill='grey',color='#708090',alpha=0.3,outlier.shape = NA,width = 0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "#708090") +
  geom_point(aes(x=subsetsize,y=pvalcount,color=subsetsize,shape=model), 
             position = position_jitter(width=0.25,height=0), size = 2,alpha=0.7) +
  scale_color_brewer(palette = 'Set1',guide = 'none')+
  scale_shape_manual(values = c(16,17,15))

ggsave("randomness_pval.pdf", plot = p, width = 3, height = 2, units = "in")


p<-ggplot(data=pvalcomb,aes(x=subsetsize,y=adjpvalcount)) +
  labs(x = "Subset Size",y = "Sig Gene Count") +
  scale_y_continuous(limits = c(750,1500), breaks=seq(from=750, to=1500, by=250)) +
  theme(plot.title=element_blank(),
        panel.background=element_rect(fill='white'),
        panel.grid.major=element_line(color='grey',linewidth =0.3),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.title=element_text(size=10),
        legend.position='top',
        legend.key.height = unit(1,'line'),
        legend.key.width  = unit(0.5,'line'),
        legend.key=element_rect(fill='transparent'),
        legend.text=element_text(size=10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.y.right = element_text(size=10),
        axis.text.y.right = element_text(size=10)) +
  geom_boxplot(fill='grey',color='#708090',alpha=0.3,outlier.shape = NA,width = 0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "#708090") +
  geom_point(aes(x=subsetsize,y=adjpvalcount,color=subsetsize,shape=model), 
             position = position_jitter(width=0.25,height=0), size = 2,alpha=0.7) +
  scale_color_brewer(palette = 'Set1',guide = 'none')+
  scale_shape_manual(values = c(16,17,15))

ggsave("randomness_adjpval.pdf", plot = p, width = 3, height = 2, units = "in")


