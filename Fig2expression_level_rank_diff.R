#load('FFPE_object.RData')
##rankdiff
TPM_rank_FF<-Lung_TPM_FF[Lung_code_TPM$Gene.ID[rowMeans(Lung_code_TPM[,-1])>10],]
TPM_rank_FF<-apply(TPM_rank_FF,2,function(x){rank(-x)})
TPM_rank_FFPE<-Lung_TPM_FFPE[Lung_code_TPM$Gene.ID[rowMeans(Lung_code_TPM[,-1])>10],]
TPM_rank_FFPE<-apply(TPM_rank_FFPE,2,function(x){rank(-x)})
TPM_rank_gap<-TPM_rank_FFPE-TPM_rank_FF
nrow(TPM_rank_gap)

dlast<-data.frame(dimnames(table(cut(abs(TPM_rank_gap[,1]), breaks = c(0,100,1000,12000)))),
                  round(as.vector(table(cut(abs(TPM_rank_gap[,1]),
                                            breaks = c(0,100,1000,12000))))/nrow(TPM_rank_gap),digits=3))
for(i in 2:15){
  dlast<-cbind(dlast,
               round(as.vector(table(cut(abs(TPM_rank_gap[,i]),
                                         breaks = c(0,100,1000,12000))))/nrow(TPM_rank_gap),digits=3))
}

colnames(dlast)<-c('Rank_difference','L7_40.5','L7_406','L7_448','L7_c0.5','L7_c06','L7_c48','L3.1','L4.1','section',
                   paste0('L',1:6))
dlast%>%gather(.,key="sample",value="proportion",-Rank_difference)%>%
  data.frame()->dlast_narrow
dlast_narrow$`Rank difference`<-factor(dlast_narrow$Rank_difference,
                                levels = c("(0,100]","(100,1e+03]","(1e+03,1.2e+04]"))
dlast_narrow$Proportion<-dlast_narrow$proportion
pdf(file = "rankdiff.pdf",height = 5,width = 9.5)
ggplot(dlast_narrow, aes( x = sample,y=Proportion,fill = `Rank difference`))+
  geom_col(position = 'stack', width = 0.6)+
  theme(legend.position='bottom',
        legend.text=element_text(size=14,colour = 'black'),
        legend.title = element_text(size=14,colour = 'black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14,angle=45,hjust=1,colour = 'black'),
        axis.title.y = element_text(size=14,colour = 'black'),
        axis.text.y = element_text(size = 14,colour = 'black'))
dev.off()