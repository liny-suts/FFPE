library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(DescTools)
library(scales)
#load('FFPE_object.RData')
##gene_attribute_information are obtained using BiomaRt and mRNAloc
Lung_nfd<-merge(Lung_deg,gene_attribute_information,by='Geneid')
Lung_nfd$down<-Lung_nfd$group=='down'
Lung_nfd$up<-Lung_nfd$group=='up'

colon_nfd<-merge(colon_deg,gene_attribute_information,by='Geneid')
colon_nfd$down<-colon_nfd$group=='down'
colon_nfd$up<-colon_nfd$group=='up'

#mark the genes
colon_nfd$exlong<-colon_nfd$transcript_length>quantile(colon_nfd$transcript_length,probs = c(0.1,0.9))[[2]]
colon_nfd$exshort<-colon_nfd$transcript_length<quantile(colon_nfd$transcript_length,probs = c(0.1,0.9))[[1]]
kidney_nfd$exlong<-kidney_nfd$transcript_length>quantile(kidney_nfd$transcript_length,probs = c(0.1,0.9))[[2]]
kidney_nfd$exshort<-kidney_nfd$transcript_length<quantile(kidney_nfd$transcript_length,probs = c(0.1,0.9))[[1]]
ovary_nfd$exlong<-ovary_nfd$transcript_length>quantile(ovary_nfd$transcript_length,probs = c(0.1,0.9))[[2]]
ovary_nfd$exshort<-ovary_nfd$transcript_length<quantile(ovary_nfd$transcript_length,probs = c(0.1,0.9))[[1]]
Lung_nfd$exlong<-Lung_nfd$transcript_length>quantile(Lung_nfd$transcript_length,probs = c(0.1,0.9))[[2]]
Lung_nfd$exshort<-Lung_nfd$transcript_length<quantile(Lung_nfd$transcript_length,probs = c(0.1,0.9))[[1]]

colon_nfd$exhigh<-colon_nfd$MFE>quantile(colon_nfd$MFE,probs = c(0.1,0.9),na.rm=T)[[2]]
colon_nfd$exlow<-colon_nfd$MFE<quantile(colon_nfd$MFE,probs = c(0.1,0.9),na.rm=T)[[1]]
kidney_nfd$exhigh<-kidney_nfd$MFE>quantile(kidney_nfd$MFE,probs = c(0.1,0.9),na.rm=T)[[2]]
kidney_nfd$exlow<-kidney_nfd$MFE<quantile(kidney_nfd$MFE,probs = c(0.1,0.9),na.rm=T)[[1]]
ovary_nfd$exhigh<-ovary_nfd$MFE>quantile(ovary_nfd$MFE,probs = c(0.1,0.9),na.rm=T)[[2]]
ovary_nfd$exlow<-ovary_nfd$MFE<quantile(ovary_nfd$MFE,probs = c(0.1,0.9),na.rm=T)[[1]]
Lung_nfd$exhigh<-Lung_nfd$MFE>quantile(Lung_nfd$MFE,probs = c(0.1,0.9),na.rm=T)[[2]]
Lung_nfd$exlow<-Lung_nfd$MFE<quantile(Lung_nfd$MFE,probs = c(0.1,0.9),na.rm=T)[[1]]

colon_nfd$exhighA<-colon_nfd$A>quantile(colon_nfd$A,probs = c(0.1,0.9),na.rm=T)[[2]]
colon_nfd$exlowA<-colon_nfd$A<quantile(colon_nfd$A,probs = c(0.1,0.9),na.rm=T)[[1]]
kidney_nfd$exhighA<-kidney_nfd$A>quantile(kidney_nfd$A,probs = c(0.1,0.9),na.rm=T)[[2]]
kidney_nfd$exlowA<-kidney_nfd$A<quantile(kidney_nfd$A,probs = c(0.1,0.9),na.rm=T)[[1]]
ovary_nfd$exhighA<-ovary_nfd$A>quantile(ovary_nfd$A,probs = c(0.1,0.9),na.rm=T)[[2]]
ovary_nfd$exlowA<-ovary_nfd$A<quantile(ovary_nfd$A,probs = c(0.1,0.9),na.rm=T)[[1]]
Lung_nfd$exhighA<-Lung_nfd$A>quantile(Lung_nfd$A,probs = c(0.1,0.9),na.rm=T)[[2]]
Lung_nfd$exlowA<-Lung_nfd$A<quantile(Lung_nfd$A,probs = c(0.1,0.9),na.rm=T)[[1]]

colon_nfd$exhighU<-colon_nfd$U>quantile(colon_nfd$U,probs = c(0.1,0.9),na.rm=T)[[2]]
colon_nfd$exlowU<-colon_nfd$U<quantile(colon_nfd$U,probs = c(0.1,0.9),na.rm=T)[[1]]
kidney_nfd$exhighU<-kidney_nfd$U>quantile(kidney_nfd$U,probs = c(0.1,0.9),na.rm=T)[[2]]
kidney_nfd$exlowU<-kidney_nfd$U<quantile(kidney_nfd$U,probs = c(0.1,0.9),na.rm=T)[[1]]
ovary_nfd$exhighU<-ovary_nfd$U>quantile(ovary_nfd$U,probs = c(0.1,0.9),na.rm=T)[[2]]
ovary_nfd$exlowU<-ovary_nfd$U<quantile(ovary_nfd$U,probs = c(0.1,0.9),na.rm=T)[[1]]
Lung_nfd$exhighU<-Lung_nfd$U>quantile(Lung_nfd$U,probs = c(0.1,0.9),na.rm=T)[[2]]
Lung_nfd$exlowU<-Lung_nfd$U<quantile(Lung_nfd$U,probs = c(0.1,0.9),na.rm=T)[[1]]

colon_nfd$exhighG<-colon_nfd$G>quantile(colon_nfd$G,probs = c(0.1,0.9),na.rm=T)[[2]]
colon_nfd$exlowG<-colon_nfd$G<quantile(colon_nfd$G,probs = c(0.1,0.9),na.rm=T)[[1]]
kidney_nfd$exhighG<-kidney_nfd$G>quantile(kidney_nfd$G,probs = c(0.1,0.9),na.rm=T)[[2]]
kidney_nfd$exlowG<-kidney_nfd$G<quantile(kidney_nfd$G,probs = c(0.1,0.9),na.rm=T)[[1]]
ovary_nfd$exhighG<-ovary_nfd$G>quantile(ovary_nfd$G,probs = c(0.1,0.9),na.rm=T)[[2]]
ovary_nfd$exlowG<-ovary_nfd$G<quantile(ovary_nfd$G,probs = c(0.1,0.9),na.rm=T)[[1]]
Lung_nfd$exhighG<-Lung_nfd$G>quantile(Lung_nfd$G,probs = c(0.1,0.9),na.rm=T)[[2]]
Lung_nfd$exlowG<-Lung_nfd$G<quantile(Lung_nfd$G,probs = c(0.1,0.9),na.rm=T)[[1]]

colon_nfd$exhighC<-colon_nfd$C>quantile(colon_nfd$C,probs = c(0.1,0.9),na.rm=T)[[2]]
colon_nfd$exlowC<-colon_nfd$C<quantile(colon_nfd$C,probs = c(0.1,0.9),na.rm=T)[[1]]
kidney_nfd$exhighC<-kidney_nfd$C>quantile(kidney_nfd$C,probs = c(0.1,0.9),na.rm=T)[[2]]
kidney_nfd$exlowC<-kidney_nfd$C<quantile(kidney_nfd$C,probs = c(0.1,0.9),na.rm=T)[[1]]
ovary_nfd$exhighC<-ovary_nfd$C>quantile(ovary_nfd$C,probs = c(0.1,0.9),na.rm=T)[[2]]
ovary_nfd$exlowC<-ovary_nfd$C<quantile(ovary_nfd$C,probs = c(0.1,0.9),na.rm=T)[[1]]
Lung_nfd$exhighC<-Lung_nfd$C>quantile(Lung_nfd$C,probs = c(0.1,0.9),na.rm=T)[[2]]
Lung_nfd$exlowC<-Lung_nfd$C<quantile(Lung_nfd$C,probs = c(0.1,0.9),na.rm=T)[[1]]
######
sig_group_violin<-function(nfd=NULL,attribute='transcript_length',log=F){
  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                      symbols = c("****", "***", "**", "*", "ns"))
a<-ggplot(data.frame(nfd), aes(x = group, y = .data[[attribute]]))+
  geom_violin(aes(fill = group),trim = FALSE)+
  theme(legend.position = 'none',
        axis.text.x = element_text(size=14,colour = 'black'),
        axis.text.y = element_text(size=14,colour = 'black'),
        axis.title.y = element_text(size=14,colour = 'black'))+
  labs(title="",x="", y = attribute)+
  geom_boxplot(width=0.3,outlier.colour = NA)+
  stat_compare_means(comparisons = list(c('down', 'ns'),c('up', 'down')),
                     symnum.args = symnum.args,)
b<-ggplot(data.frame(nfd), aes(x = group, y = log2(.data[[attribute]])))+
  geom_violin(aes(fill = group),trim = FALSE)+
  theme(legend.position = 'none',
        axis.text.x = element_text(size=14,colour = 'black'),
        axis.text.y = element_text(size=14,colour = 'black'),
        axis.title.y = element_text(size=14,colour = 'black'))+
  labs(title="",x="", y = paste0("log2 ",attribute))+
  geom_boxplot(width=0.3,outlier.colour = NA)+
  stat_compare_means(comparisons = list(c('down', 'ns'),c('up', 'down')),
                     symnum.args = symnum.args)
if(log){return(b)}
else{return(a)}
}
###fisher exact test visualization
nfanalyze<-function(nfd=NULL,sig='up',ydown=-7,yup=5){
dplyr::count(nfd,.data[[sig]],exlong)%>%tidyr::spread(key=exlong,value = n)%>%data.frame()->tmp
tmp[is.na(tmp)]<-1
fisher.test(tmp[,c(-1)])->Plong
dplyr::count(nfd,.data[[sig]],exshort)%>%tidyr::spread(key=exshort,value = n)%>%data.frame()->tmp
tmp[is.na(tmp)]<-1
fisher.test(tmp[,c(-1)])->Pshort
nfd%>%filter(.,!is.na(MFE))%>%dplyr::count(.,.data[[sig]],exhigh)%>%tidyr::spread(key=exhigh,value = n)%>%data.frame()->tmp
tmp[is.na(tmp)]<-1
fisher.test(tmp[,c(-1)])->Phigh
nfd%>%filter(.,!is.na(MFE))%>%dplyr::count(.,.data[[sig]],exlow)%>%tidyr::spread(key=exlow,value = n)%>%data.frame()->tmp
tmp[is.na(tmp)]<-1
fisher.test(tmp[,c(-1)])->Plow
dplyr::count(nfd,.data[[sig]],exhighG)%>%tidyr::spread(key=exhighG,value = n)%>%data.frame()->tmp
tmp[is.na(tmp)]<-1
fisher.test(tmp[,c(-1)])->PhighG
dplyr::count(nfd,.data[[sig]],exlowG)%>%tidyr::spread(key=exlowG,value = n)%>%data.frame()->tmp
tmp[is.na(tmp)]<-1
fisher.test(tmp[,c(-1)])->PlowG
dplyr::count(nfd,.data[[sig]],exhighC)%>%tidyr::spread(key=exhighC,value = n)%>%data.frame()->tmp
tmp[is.na(tmp)]<-1
fisher.test(tmp[,c(-1)])->PhighC
dplyr::count(nfd,.data[[sig]],exlowC)%>%tidyr::spread(key=exlowC,value = n)%>%data.frame()->tmp
tmp[is.na(tmp)]<-1
fisher.test(tmp[,c(-1)])->PlowC
dplyr::count(nfd,.data[[sig]],exhighA)%>%tidyr::spread(key=exhighA,value = n)%>%data.frame()->tmp
tmp[is.na(tmp)]<-1
fisher.test(tmp[,c(-1)])->PhighA
dplyr::count(nfd,.data[[sig]],exlowA)%>%tidyr::spread(key=exlowA,value = n)%>%data.frame()->tmp
tmp[is.na(tmp)]<-1
fisher.test(tmp[,c(-1)])->PlowA
dplyr::count(nfd,.data[[sig]],exhighU)%>%tidyr::spread(key=exhighU,value = n)%>%data.frame()->tmp
tmp[is.na(tmp)]<-1
fisher.test(tmp[,c(-1)])->PhighU
dplyr::count(nfd,.data[[sig]],exlowU)%>%tidyr::spread(key=exlowU,value = n)%>%data.frame()->tmp
tmp[is.na(tmp)]<-1
fisher.test(tmp[,c(-1)])->PlowU
dplyr::count(nfd,.data[[sig]],cytoplasm)%>%tidyr::spread(key=cytoplasm,value = n)%>%data.frame()->tmp
fisher.test(tmp[,c(-1)])->Pcytoplasm
dplyr::count(nfd,.data[[sig]],ER)%>%tidyr::spread(key=ER,value = n)%>%data.frame()->tmp
fisher.test(tmp[,c(-1)])->PER
dplyr::count(nfd,.data[[sig]],nucleus)%>%tidyr::spread(key=nucleus,value = n)%>%data.frame()->tmp
fisher.test(tmp[,c(-1)])->Pnucleus
dplyr::count(nfd,.data[[sig]],extracellular)%>%tidyr::spread(key=extracellular,value = n)%>%data.frame()->tmp
fisher.test(tmp[,c(-1)])->Pextracellular

downf_e_tres<-data.frame('length>D9'=c(Plong$estimate,Plong$conf.int[1]+0.0075,Plong$conf.int[2],signif(Plong$p.value,2)),
                         'length<D1'=c(Pshort$estimate,Pshort$conf.int[1]+0.0075,Pshort$conf.int[2],signif(Pshort$p.value,2)),
                         'MFE>D9'=c(Phigh$estimate,Phigh$conf.int[1]+0.0075,Phigh$conf.int[2],signif(Phigh$p.value,2)),
                         'MFE<D1'=c(Plow$estimate,Plow$conf.int[1]+0.0075,Plow$conf.int[2],signif(Plow$p.value,2)),
                         'G%>D9'=c(PhighG$estimate,PhighG$conf.int[1]+0.0075,PhighG$conf.int[2],signif(PhighG$p.value,2)),
                         'G%<D1'=c(PlowG$estimate,PlowG$conf.int[1]+0.0075,PlowG$conf.int[2],signif(PlowG$p.value,2)),
                         'C%>D9'=c(PhighC$estimate,PhighC$conf.int[1]+0.0075,PhighC$conf.int[2],signif(PhighC$p.value,2)),
                         'C%<D1'=c(PlowC$estimate,PlowC$conf.int[1]+0.0075,PlowC$conf.int[2],signif(PlowC$p.value,2)),
                         'A%>D9'=c(PhighA$estimate,PhighA$conf.int[1]+0.0075,PhighA$conf.int[2],signif(PhighA$p.value,2)),
                         'A%<D1'=c(PlowA$estimate,PlowA$conf.int[1]+0.0075,PlowA$conf.int[2],signif(PlowA$p.value,2)),
                         'U%>D9'=c(PhighU$estimate,PhighU$conf.int[1]+0.0075,PhighU$conf.int[2],signif(PhighU$p.value,2)),
                         'U%<D1'=c(PlowU$estimate,PlowU$conf.int[1]+0.0075,PlowU$conf.int[2],signif(PlowU$p.value,2)),
                         'cytoplasm'=c(Pcytoplasm$estimate,Pcytoplasm$conf.int[1],Pcytoplasm$conf.int[2],signif(Pcytoplasm$p.value,2)),
                         'ER'=c(PER$estimate,PER$conf.int[1],PER$conf.int[2],signif(PER$p.value,2)),
                         'nucleus'=c(Pnucleus$estimate,Pnucleus$conf.int[1],Pnucleus$conf.int[2],signif(Pnucleus$p.value,2)),
                         'extracellular'=c(Pextracellular$estimate,Pextracellular$conf.int[1],Pextracellular$conf.int[2],signif(Pextracellular$p.value,2)))
downf_e_tres<-data.frame(t(downf_e_tres))
downf_e_tres$color<-ifelse(abs(log2(downf_e_tres[,1]))>1&downf_e_tres[,4]<0.005,'red','black')
colnames(downf_e_tres)[1:4]<-c('odds_ratio','Down_CI','Up_CI','pvalue')
downf_e_tres$geneset<-factor(c('length>D9','length<D1','MFE>D9','MFE<D1','G%>D9','G%<D1','C%>D9','C%<D1','A%>D9','A%<D1','U%>D9','U%<D1','cytoplasm','ER','nucleus','extracellular'),
                             levels=c('length>D9','length<D1','MFE>D9','MFE<D1','G%>D9','G%<D1','C%>D9','C%<D1','A%>D9','A%<D1','U%>D9','U%<D1','cytoplasm','ER','nucleus','extracellular'))
group <- c(red = "red",black = "grey")
downf_e_tres$`log2(odds ratio)`<-log2(downf_e_tres$odds_ratio)
downf_e_tres$Geneset<-downf_e_tres$geneset
a<-ggplot(downf_e_tres,aes(x=Geneset,color=color))+
  scale_color_manual(values = group)+
  geom_errorbar(aes(x=Geneset,ymin=log2(Down_CI),ymax=log2(Up_CI),color=color),width=.3)+
  geom_point(aes(x=Geneset,y=`log2(odds ratio)`,color=color),size=3)+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.position = 'none',axis.text.x = element_text(size=16,color = 'black'),
        axis.text.y = element_text(size=16,color = 'black'),axis.title.x = element_text(size=15,color = 'black'),
        axis.title.y = element_text(size=15,color = 'black'))+geom_hline(yintercept=c(0,-1,1), size=c(1,0.5,0.5),linetype="dotted")+
  ylim(ydown,yup)+
  coord_flip()
tmp<-data.frame(apply(downf_e_tres[,1:3],2,function(x){return(format(round(x,2),digits=1,nsmall=2))}))
newtmp<-data.frame('Odds_ratio(CI:95%)'=paste0(Rev(tmp$odds_ratio),'  [',Rev(tmp$Down_CI),',',Rev(tmp$Up_CI),']'),'Pvalue'=Rev(downf_e_tres$pvalue))
colnames(newtmp)<-c('Odds_ratio(CI:95%)','p_value')
newtmp$a<-newtmp$p_value/10^floor(log10(newtmp$p_value))
newtmp$b<-floor(log10(newtmp$p_value))

b<-ggplot(newtmp) +
  annotate("text", x = 0, y = 16:1, label = newtmp$`Odds_ratio(CI:95%)`,size=5)+
  annotate("text", x = 0.002, y = 16, label = bquote(.(newtmp$a[1])%*%10^.(newtmp$b[1])),size=5)+
  annotate("text", x = 0.002, y = 15, label = bquote(.(newtmp$a[2])%*%10^.(newtmp$b[2])),size=5)+
  annotate("text", x = 0.002, y = 14, label = bquote(.(newtmp$a[3])%*%10^.(newtmp$b[3])),size=5)+
  annotate("text", x = 0.002, y = 13, label = bquote(.(newtmp$a[4])%*%10^.(newtmp$b[4])),size=5)+
  annotate("text", x = 0.002, y = 12, label = bquote(.(newtmp$a[5])%*%10^.(newtmp$b[5])),size=5)+
  annotate("text", x = 0.002, y = 11, label = bquote(.(newtmp$a[6])%*%10^.(newtmp$b[6])),size=5)+
  annotate("text", x = 0.002, y = 10, label = bquote(.(newtmp$a[7])%*%10^.(newtmp$b[7])),size=5)+
  annotate("text", x = 0.002, y = 9, label = bquote(.(newtmp$a[8])%*%10^.(newtmp$b[8])),size=5)+
  annotate("text", x = 0.002, y = 8, label = bquote(.(newtmp$a[9])%*%10^.(newtmp$b[9])),size=5)+
  annotate("text", x = 0.002, y = 7, label = bquote(.(newtmp$a[10])%*%10^.(newtmp$b[10])),size=5)+
  annotate("text", x = 0.002, y = 6, label = bquote(.(newtmp$a[11])%*%10^.(newtmp$b[11])),size=5)+
  annotate("text", x = 0.002, y = 5, label = bquote(.(newtmp$a[12])%*%10^.(newtmp$b[12])),size=5)+
  annotate("text", x = 0.002, y = 4, label = bquote(.(newtmp$a[13])%*%10^.(newtmp$b[13])),size=5)+
  annotate("text", x = 0.002, y = 3, label = bquote(.(newtmp$a[14])%*%10^.(newtmp$b[14])),size=5)+
  annotate("text", x = 0.002, y = 2, label = bquote(.(newtmp$a[15])%*%10^.(newtmp$b[15])),size=5)+
  annotate("text", x = 0.002, y = 1, label = bquote(.(newtmp$a[16])%*%10^.(newtmp$b[16])),size=5)+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  ylab('')+xlab('')+xlim(-0.0009,0.0028)
return(a|b)
}
pdf("up_nfd_colon.pdf",width = 7.5,height = 4.8)
nfanalyze(nfd = colon_nfd,sig = 'up',ydown = -6)
dev.off()
pdf("up_nfd_kidney.pdf",width = 7.5,height = 4.8)
nfanalyze(nfd = kidney_nfd,sig = 'up',ydown = -6)
dev.off()
pdf("up_nfd_ovary.pdf",width = 7.5,height = 4.8)
nfanalyze(nfd = ovary_nfd,sig = 'up',ydown = -6)
dev.off()
pdf("up_nfd_lung.pdf",width = 7.5,height = 4.8)
nfanalyze(nfd = Lung_nfd,sig = 'up',ydown = -6)
dev.off()

pdf("down_nfd_colon.pdf",width = 8.4,height = 4.8)
nfanalyze(nfd = colon_nfd,sig = 'down',ydown = -8)
dev.off()
pdf("down_nfd_kidney.pdf",width = 7.5,height = 4.8)
nfanalyze(nfd = kidney_nfd,sig = 'down',ydown = -8)
dev.off()
pdf("down_nfd_ovary.pdf",width = 8.4,height = 4.8)
nfanalyze(nfd = ovary_nfd,sig = 'down',ydown = -8)
dev.off()
pdf("down_nfd_lung.pdf",width = 7.5,height = 4.8)
nfanalyze(nfd = Lung_nfd,sig = 'down',ydown = -8)
dev.off()

pdf("transcript_length.pdf",width = 5,height = 4.1)
sig_group_violin(nfd = colon_nfd,attribute = 'transcript_length',log=T)
dev.off()

pdf("MFE.pdf",width = 5,height = 4.1)
sig_group_violin(nfd = colon_nfd,attribute = 'MFE')
dev.off()

pdf("Lungtranscript_length.pdf",width = 4.5,height = 3.5)
sig_group_violin(nfd = Lung_nfd,attribute = 'transcript_length',log=T)+
  ylab('log2(transcript length)')
dev.off()

pdf("LungMFE.pdf",width = 5,height = 3.5)
sig_group_violin(nfd = Lung_nfd,attribute = 'MFE')+ylim(-7500,2000)
dev.off()

pdf("kidneytranscript_length.pdf",width = 5,height = 4.1)
sig_group_violin(nfd = kidney_nfd,attribute = 'transcript_length',log=T)
dev.off()

pdf("kidneyMFE.pdf",width = 5,height = 4.1)
sig_group_violin(nfd = kidney_nfd,attribute = 'MFE')
dev.off()

pdf("ovarytranscript_length.pdf",width = 5,height = 4.1)
sig_group_violin(nfd = ovary_nfd,attribute = 'transcript_length',log=T)
dev.off()

pdf("ovaryMFE.pdf",width = 5,height = 4.1)
sig_group_violin(nfd = ovary_nfd,attribute = 'MFE')
dev.off()

