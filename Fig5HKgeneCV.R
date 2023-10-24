library(ggplot2)
library(ggpubr)
library(tidyverse)
#load('FFPE_object.RData')
##Housekeeping gene information was obtain from HRT atlas.
cv_1<-Lung_TPM_FF[FHK$Gene.ID,c(7,8,10:15)]%>%
  apply(.,1,function(x){sd(x)/mean(x+0.1)})
cv_2<-Lung_TPM_FFPE[FHK$Gene.ID,c(7,8,10:15)]%>%
  apply(.,1,function(x){sd(x)/mean(x+0.1)})

BPVpaired_TPM[FHK$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='OVARY'&BPVpaired_meta$PROTOCOL=='CI_1h']]%>%
  apply(.,1,function(x){sd(x)/mean(x+0.1)})->cv_O1h
BPVpaired_TPM[FHK$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='OVARY'&BPVpaired_meta$PROTOCOL=='CI_2h']]%>%
  apply(.,1,function(x){sd(x)/mean(x+0.1)})->cv_O2h
BPVpaired_TPM[FHK$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='OVARY'&BPVpaired_meta$PROTOCOL=='CI_3h']]%>%
  apply(.,1,function(x){sd(x)/mean(x+0.1)})->cv_O3h
BPVpaired_TPM[FHK$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='OVARY'&BPVpaired_meta$PROTOCOL=='CI_12h']]%>%
  apply(.,1,function(x){sd(x)/mean(x+0.1)})->cv_O12h
BPVpaired_TPM[FHK$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='OVARY'&BPVpaired_meta$PROTOCOL=='Frozen']]%>%
  apply(.,1,function(x){sd(x)/mean(x+0.1)})->cv_OFF

BPVpaired_TPM[FHK$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='COLON'&BPVpaired_meta$PROTOCOL=='CI_1h']]%>%
  apply(.,1,function(x){sd(x)/mean(x+0.1)})->cv_C1h
BPVpaired_TPM[FHK$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='COLON'&BPVpaired_meta$PROTOCOL=='CI_2h']]%>%
  apply(.,1,function(x){sd(x)/mean(x+0.1)})->cv_C2h
BPVpaired_TPM[FHK$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='COLON'&BPVpaired_meta$PROTOCOL=='CI_3h']]%>%
  apply(.,1,function(x){sd(x)/mean(x+0.1)})->cv_C3h
BPVpaired_TPM[FHK$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='COLON'&BPVpaired_meta$PROTOCOL=='CI_12h']]%>%
  apply(.,1,function(x){sd(x)/mean(x+0.1)})->cv_C12h
BPVpaired_TPM[FHK$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='COLON'&BPVpaired_meta$PROTOCOL=='Frozen']]%>%
  apply(.,1,function(x){sd(x)/mean(x+0.1)})->cv_CFF

BPVpaired_TPM[FHK$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='KIDNEY'&BPVpaired_meta$PROTOCOL=='CI_1h']]%>%
  apply(.,1,function(x){sd(x)/mean(x+0.1)})->cv_K1h
BPVpaired_TPM[FHK$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='KIDNEY'&BPVpaired_meta$PROTOCOL=='CI_2h']]%>%
  apply(.,1,function(x){sd(x)/mean(x+0.1)})->cv_K2h
BPVpaired_TPM[FHK$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='KIDNEY'&BPVpaired_meta$PROTOCOL=='CI_3h']]%>%
  apply(.,1,function(x){sd(x)/mean(x+0.1)})->cv_K3h
BPVpaired_TPM[FHK$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='KIDNEY'&BPVpaired_meta$PROTOCOL=='CI_12h']]%>%
  apply(.,1,function(x){sd(x)/mean(x+0.1)})->cv_K12h
BPVpaired_TPM[FHK$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='KIDNEY'&BPVpaired_meta$PROTOCOL=='Frozen']]%>%
  apply(.,1,function(x){sd(x)/mean(x+0.1)})->cv_KFF

cv_frame<-data.frame("Lung_FF"=cv_1,"Lung_FFPE"=cv_2,
           "Colon_FF"=cv_CFF,"Colon_1h"=cv_C1h,"Colon_2h"=cv_C2h,"Colon_3h"=cv_C3h,"Colon_12h"=cv_C12h,
           "Kidney_FF"=cv_KFF,"Kidney_1h"=cv_K1h,"Kidney_2h"=cv_K2h,"Kidney_3h"=cv_K3h,"Kidney_12h"=cv_K12h,
           "Ovary_FF"=cv_OFF,"Ovary_1h"=cv_O1h,"Ovary_2h"=cv_O2h,"Ovary_3h"=cv_O3h,"Ovary_12h"=cv_O12h)

tmp<-cv_frame%>%gather(.,key = "group",value = "CV")
tmp$group<-factor(tmp$group,levels = unique(tmp$group))
pdf("/home/liny/R/labworkbase/FFPE_final/FFPE_final_plot/plot/HKgene/CVplot.pdf",width = 7,height = 5)
ggplot(tmp, aes(x = group, y = CV)) +
  geom_violin(trim = FALSE, scale = "width", fill = "lightblue", alpha = 0.5) +
  geom_boxplot(width = 0.3, fill = "white", color = "black", alpha = 0.7) +
  theme_classic() +
  labs(title = "Coefficient of Variation of HK Genes",
       x = "", y = "Coefficient of Variation") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(color = 'black',angle = 45,hjust = 1),
        axis.text.y = element_text(color = 'black'))
dev.off()

##check gene whether express
cv_frame<-do.call(cbind, lapply(lapply(cv_frame, unlist), `length<-`, max(lengths(cv_frame))))
tmp<-cv_frame%>%data.frame()%>%gather(.,key = "group",value = "CV")
tmp<-tmp[!is.na(tmp$CV),]
tmp$group<-factor(tmp$group,levels = unique(tmp$group))
pdf("GAPDH.pdf",width = 5,height = 5)
ggplot(tmp, aes(x = group, y = log2(CV),fill = group)) +
  geom_boxplot(width = 0.8, color = "black", alpha = 1) +
  theme_classic() +
  labs(title = "GAPDH",
       x = "", y = "log2(TPM)") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10,color="black"),
        axis.text.x = element_text(angle = 45,hjust = 1,color = "black"),)+
  guides(fill=F)+
  stat_compare_means(comparisons = list(c("Lung_FF","Lung_FFPE"),c("Colon_FF","Colon_FFPE"),
                                        c("Kidney_FF","Kidney_FFPE"),c("Ovary_FF","Ovary_FFPE")),paired = T)
dev.off()

#check gene whether express
tmp<-c(
  sum(apply(BPV_TPM_FF[FHK1$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='COLON'&BPVpaired_meta$FF=='FFPE']],1,min)>10),
  sum(apply(BPV_TPM_FFPE[FHK1$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='COLON'&BPVpaired_meta$FF=='FFPE']],1,min)>10),
  sum(apply(BPV_TPM_FF[FHK1$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='KIDNEY'&BPVpaired_meta$FF=='FFPE']],1,min)>10),
  sum(apply(BPV_TPM_FFPE[FHK1$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='KIDNEY'&BPVpaired_meta$FF=='FFPE']],1,min)>10),
  sum(apply(Lung_TPM_FF[FHK1$Gene.ID,],1,min)>10),
  sum(apply(Lung_TPM_FFPE[FHK1$Gene.ID,],1,min)>10),
  sum(apply(BPV_TPM_FF[FHK1$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='OVARY'&BPVpaired_meta$FF=='FFPE']],1,min)>10),
  sum(apply(BPV_TPM_FFPE[FHK1$Gene.ID,BPVpaired_meta$Run[BPVpaired_meta$Body.Site=='OVARY'&BPVpaired_meta$FF=='FFPE']],1,min)>10)
)
tmp1<-data.frame('ngenes'=tmp,
                 'type'=rep(c('FF','FFPE'),4),
                 'site'=c('Colon','Colon','Kidney','Kidney','Lung','Lung','Ovary','Ovary'))
pdf("HKexpresssummary.pdf",width = 7,height = 5)
ggplot(tmp1,aes(x=type,y=ngenes,fill=type))+
  geom_col(width=0.6)+
  geom_text(aes(label=ngenes, y=ngenes+1), position=position_dodge(0.9), vjust=0)+
  facet_wrap(~site,nrow = 1)+
  guides(fill=F)+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text.y = element_text(size=12,color="black"),
        axis.text.x = element_text(size=12,color="black"),
        axis.title.y = element_text(size=12,color="black"),
        axis.title.x = element_text(size=12,color="black"),)+
  labs(title='HK genes (TPM>10)',x='',y='Number of genes')+
  theme(text=element_text(size=12))+
  ylim(0,2200)
dev.off()