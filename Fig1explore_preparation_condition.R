library(tidyverse)
library(ggplot2)
library(ggpubr)
library(table1)
library(gt)
library(DescTools)
#load('FFPE_object.RData')
##fixation time and temperature
fixation_RNA_quality <- read.csv("data_in_supplementary1")
fixation_RNA_quality$time<-factor(fixation_RNA_quality$time,levels = unique(fixation_RNA_quality$time))
fixation_RNA_quality$temperature<-factor(fixation_RNA_quality$temperature,levels = unique(fixation_RNA_quality$temperature))

pdf("DV200_fixation_time.pdf",height = 5,width = 4.5)
ggplot(fixation_RNA_quality,aes(x=time,y=DV200,fill=time))+
  geom_boxplot()+
  geom_jitter(width = 0.05)+
  theme_bw()+theme(panel.grid.minor=element_blank())+
  theme(axis.text.y = element_text(size=14,color="black"),
        axis.text.x = element_text(size=14,color="black"),
        axis.title.y = element_text(size=14,color="black"),
        axis.title.x = element_text(size=14,color="black"),
        panel.grid.major.x = element_blank())+
  guides(fill=F)+
  xlab("Fixation time")+
  stat_compare_means(comparisons = list(c('12h','24h'),c('12h','48h'),c('24h','48h')),paired = T)
dev.off()

pdf("DV800_fixation_time.pdf",height = 5,width = 4.5)
ggplot(fixation_RNA_quality,aes(x=time,y=DV800,fill=time))+
  geom_boxplot()+
  geom_jitter(width = 0.05)+
  theme_bw()+theme(panel.grid.minor=element_blank())+
  theme(axis.text.y = element_text(size=14,color="black"),
        axis.text.x = element_text(size=14,color="black"),
        axis.title.y = element_text(size=14,color="black"),
        axis.title.x = element_text(size=14,color="black"),
        panel.grid.major.x = element_blank())+
  guides(fill=F)+
  xlab("Fixation time")+
  stat_compare_means(comparisons = list(c('12h','24h'),c('12h','48h'),c('24h','48h')),paired = T)
dev.off()

pdf("DV200_fixation_temperature.pdf",height = 5,width = 3.5)
ggplot(fixation_RNA_quality[,],aes(x=temperature,y=DV200,fill=temperature))+
  geom_boxplot()+
  geom_jitter(width = 0.05)+
  theme_bw()+theme(panel.grid.minor=element_blank())+
  theme(axis.text.y = element_text(size=14,color="black"),
        axis.text.x = element_text(size=14,color="black"),
        axis.title.y = element_text(size=14,color="black"),
        axis.title.x = element_text(size=14,color="black"),
        panel.grid.major.x = element_blank())+
  guides(fill=F)+
  xlab("Fixation temperature")+
  stat_compare_means(paired = T,label.x = 1.4)
dev.off()

pdf("DV800_fixation_temperature.pdf",height = 5,width = 3.5)
ggplot(fixation_RNA_quality[,],aes(x=temperature,y=DV800,fill=temperature))+
  geom_boxplot()+
  geom_jitter(width = 0.05)+
  theme_bw()+theme(panel.grid.minor=element_blank())+
  theme(axis.text.y = element_text(size=14,color="black"),
        axis.text.x = element_text(size=14,color="black"),
        axis.title.y = element_text(size=14,color="black"),
        axis.title.x = element_text(size=14,color="black"),
        panel.grid.major.x = element_blank())+
  guides(fill=F)+
  xlab("Fixation temperature")+
  stat_compare_means(paired = T,label.x = 1.4)
dev.off()

## ishemia time
CIT_time <- read.csv("CIT_time.csv")
CIT_time$CIT<-factor(CIT_time$CIT,levels = unique(CIT_time$CIT))

pdf("DV200_CIT.pdf",height = 5,width = 4.5)
ggplot(CIT_time,aes(x=CIT,y=DV200,fill=CIT))+
  geom_boxplot()+
  geom_jitter(width = 0.05)+
  theme_bw()+theme(panel.grid.minor=element_blank())+
  theme(axis.text.y = element_text(size=14,color="black"),
        axis.text.x = element_text(size=14,color="black"),
        axis.title.y = element_text(size=14,color="black"),
        axis.title.x = element_text(size=14,color="black"),
        panel.grid.major.x = element_blank())+
  guides(fill=F)+
  xlab("Cold ischemia time")+
  stat_compare_means(paired = T,size=5)#comparisons = list(c('0.5h','3h'),c('0.5h','6h'),c('0.5h','12h'))
dev.off()

##CCCanalysis

get_CCC<-function(x,y){
  corres<-CCC(as.matrix(x[,1]),as.matrix(y[,1]))
  corres<-corres$rho.c
  for (i in 2:ncol(x)) {
    CCC(as.matrix(x[,i]),as.matrix(y[,i]))->corres1
    corres1<-corres1$rho.c
    corres<-rbind(corres,corres1)
  }
  rownames(corres)<-c(1:ncol(x))
  return(corres)
}
P_CCCres<-get_CCC(Lung_TPM_FF,Lung_TPM_FFPE)
P_CCCres<-data.frame(P_CCCres,row.names = gsub("FFPE\\.","",colnames(Lung_TPM_FF)))
P_CCCres<-data.frame(apply(P_CCCres,2,function(x){return(format(x,digit=1,nsmall=2))}))
P_ranked_gene<-apply(cbind(Lung_TPM_FF,Lung_TPM_FFPE),2,function(x){return(rownames(Lung_TPM_FF)[order(x,decreasing = T)])})
top_overlap<-function(x,y,n=10){
  res<-c(length(intersect(as.matrix(x[1:n,1]),as.matrix(y[1:n,1]))))
  for (i in 2:ncol(x)) {
    res<-append(res,length(intersect(as.matrix(x[1:n,i]),as.matrix(y[1:n,i]))))
  }
  return(res)
}
newtmp<-tibble('10'=top_overlap(P_ranked_gene[,1:15],P_ranked_gene[,16:30]),
               '100'=top_overlap(P_ranked_gene[,1:15],P_ranked_gene[,16:30],n=100),
               '1000'=top_overlap(P_ranked_gene[,1:15],P_ranked_gene[,16:30],n=1000),
               'rho.c'=P_CCCres$est,
               '(95%CI)'=paste0('[',P_CCCres$lwr.ci,',',P_CCCres$upr.ci,']'),
               row.names=rownames(P_CCCres))
write.csv(newtmp,"CCCanalysis.csv")
rm(newtmp)

#CIT condition col plot
(CIT_condtion<-newtmp[c(1:6),c(4,6)])
CIT_condtion$rho.c<-as.numeric(CIT_condtion$rho.c)
CIT_condtion$condition<-factor(c('4°C_0.5h','4°C_6h','4°C_48h','25°C_0.5h','25°C_6h','25°C_48h'),
                               levels = c('4°C_0.5h','4°C_6h','4°C_48h','25°C_0.5h','25°C_6h','25°C_48h'))
CIT_condtion$tmp<-c(rep('4',3),rep('25',3))
pdf("CI condition.pdf",height = 4.5,width = 5)
ggplot(CIT_condtion,aes(x=condition,y=rho.c,fill=tmp))+
  geom_col()+
  theme_bw()+theme(panel.grid.minor=element_blank())+
  theme(axis.text.y = element_text(size=14,color="black"),
        axis.text.x = element_text(size=14,color="black",angle = 30,hjust=1),
        axis.title.y = element_text(size=14,color="black"),
        axis.title.x = element_text(size=14,color="black"),
        panel.grid.major.x = element_blank())+
  guides(fill=F)+
  xlab("Cold ischemic condition")+
  ylab("CCC")
dev.off()
sampling$`Sampling method`<-sampling$type
pdf("sampling_method.pdf",height = 3.5,width = 5)
ggplot(sampling,aes(x=CIT,y=DV200,fill=`Sampling method`))+
  geom_line(mapping = aes(color=`Sampling method`),size=2)+
  theme_bw()+theme(panel.grid.minor=element_blank())+
  theme(axis.text.y = element_text(size=14,color="black"),
        axis.text.x = element_text(size=14,color="black"),
        axis.title.y = element_text(size=14,color="black"),
        axis.title.x = element_text(size=14,color="black"),
        legend.position = 'bottom',
        legend.title = element_text(size=12,color = 'black'),
        legend.text = element_text(size=12,color = 'black'))
dev.off()

pdf("sampling_method_col.pdf",height = 4,width = 7)
ggplot(sampling,aes(x=CIT,y=DV200,fill=`Sampling method`))+
  geom_bar(stat = 'identity')+
  facet_grid(.~`Sampling method`,)+
  theme_bw()+theme(panel.grid.minor=element_blank())+
  theme(axis.text.y = element_text(size=14,color="black"),
        axis.text.x = element_text(size=14,color="black"),
        axis.title.y = element_text(size=14,color="black"),
        axis.title.x = element_text(size=14,color="black"),
        legend.position = 'none',
        legend.title = element_text(size=12,color = 'black'),
        legend.text = element_text(size=12,color = 'black'))
dev.off()
