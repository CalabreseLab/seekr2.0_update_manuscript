# barplot represent venn diagram
# set working directory to where the files are located


k6<-read.csv('v43_v_XIST_pval_6.csv',header=T)
k5<-read.csv('v43_v_XIST_pval_5.csv',header=T)
k4<-read.csv('v43_v_XIST_pval_4.csv',header=T)

k6<-k6$X[which(k6$ENST00000429829.6_XIST<0.05)] # 2769
k5<-k5$X[which(k5$ENST00000429829.6_XIST<0.05)] # 2817
k4<-k4$X[which(k4$ENST00000429829.6_XIST<0.05)] # 2570

# needs to remove XIST itself from the list
k6<-k6[!grepl('XIST',k6)] # 2768
k5<-k5[!grepl('XIST',k5)] # 2816
k4<-k4[!grepl('XIST',k4)] # 2569


df<-data.frame(group=c('K-456','K-56','K-46','K-45',
                       'K-6','K-5','K-4'),
               count=rep(0,times=7))

df$count[which(df$group=='K-456')]<-length(intersect(intersect(k6,k5),k4))

df$count[which(df$group=='K-56')]<-length(intersect(k6,k5))-df$count[which(df$group=='K-456')]

df$count[which(df$group=='K-46')]<-length(intersect(k6,k4))-df$count[which(df$group=='K-456')]

df$count[which(df$group=='K-45')]<-length(intersect(k5,k4))-df$count[which(df$group=='K-456')]

df$count[which(df$group=='K-6')]<-length(k6)-df$count[which(df$group=='K-456')]-
  df$count[which(df$group=='K-56')]-df$count[which(df$group=='K-46')]

df$count[which(df$group=='K-5')]<-length(k5)-df$count[which(df$group=='K-456')]-
  df$count[which(df$group=='K-56')]-df$count[which(df$group=='K-45')]

df$count[which(df$group=='K-4')]<-length(k4)-df$count[which(df$group=='K-456')]-
  df$count[which(df$group=='K-46')]-df$count[which(df$group=='K-45')]


df$group<-factor(df$group,levels=c('K-456','K-56','K-46','K-45',
                                   'K-6','K-5','K-4'))

k456_uni<-sum(df$count)

df$percentage<-df$count*100/k456_uni

library(ggplot2)


p<-ggplot(df, aes(x=group, y=percentage))+
  scale_x_discrete('Group')+
  scale_y_continuous('% of total in Group')+
  theme(plot.title=element_blank(),
        panel.background=element_rect(fill='white'),
        panel.grid.major=element_line(color='grey',linewidth =0.3),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.text.x=element_text(size=24,face='bold',color='black',angle=45,vjust = 1, hjust=1),
        axis.text.y=element_text(size=24,face='bold',color='black'),
        axis.title.x=element_text(size=28,face='bold'),
        axis.title.y=element_text(size=28,face='bold',angle=90))+
  geom_bar(stat="identity",fill='darkgrey',color='black',linewidth=1)

ggsave("barplot_pval_k456.pdf", plot = p, width = 12, height = 8, units = "in")

