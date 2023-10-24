library(package = "org.Hs.eg.db")
library(package = "clusterProfiler")
library(package = "GSEABase")
library(package = GSVA)
library(ggplot2)
library(ggpubr)
library(stringr)
library(biomaRt)
#load('FFPE_object.RData')
###GSEA
filter=rownames(ovary_deg)[ovary_deg$baseMean>1]
genelist<-ovary_deg[filter,2]
names(genelist)<-rownames(ovary_deg[filter,])
genelist<-sort(genelist,decreasing = T)
GsGO<-gseGO(
  genelist,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = T,
  by = "fgsea"
)
GsGOres<-GsGO@result[order(GsGO@result$qvalue,decreasing = F),]
GOBPdown<-GsGOres[GsGOres$NES<0,]
GOBPup<-GsGOres[GsGOres$NES>0,]
ovaryGSEAtop10up<-GOBPup[1:10,]
ovaryGSEAtop10down<-GOBPdown[1:10,]

filter=rownames(colon_deg)[colon_deg$baseMean>1]
genelist<-colon_deg[filter,2]
names(genelist)<-rownames(colon_deg[filter,])
genelist<-sort(genelist,decreasing = T)
GsGO<-gseGO(
  genelist,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = T,
  by = "fgsea"
)
GsGOres<-GsGO@result[order(GsGO@result$qvalue,decreasing = F),]
GOBPdown<-GsGOres[GsGOres$NES<0,]
GOBPup<-GsGOres[GsGOres$NES>0,]
colonGSEAtop10up<-GOBPup[1:10,]
colonGSEAtop10down<-GOBPdown[1:10,]

filter=rownames(kidney_deg)[kidney_deg$baseMean>1]
genelist<-kidney_deg[filter,2]
names(genelist)<-rownames(kidney_deg[filter,])
genelist<-sort(genelist,decreasing = T)
GsGO<-gseGO(
  genelist,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = T,
  by = "fgsea"
)
GsGOres<-GsGO@result[order(GsGO@result$qvalue,decreasing = F),]
GOBPdown<-GsGOres[GsGOres$NES<0,]
GOBPup<-GsGOres[GsGOres$NES>0,]
kidneyGSEAtop10up<-GOBPup[1:10,]
kidneyGSEAtop10down<-GOBPdown[1:10,]

filter=rownames(Lung_deg)[Lung_deg$baseMean>1]
genelist<-Lung_deg[filter,2]
names(genelist)<-rownames(Lung_deg[filter,])
genelist<-sort(genelist,decreasing = T)
GsGO<-gseGO(
  genelist,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = T,
  by = "fgsea"
)
GsGOres<-GsGO@result[order(GsGO@result$qvalue,decreasing = F),]
GOBPdown<-GsGOres[GsGOres$NES<0,]
GOBPup<-GsGOres[GsGOres$NES>0,]
LungGSEAtop10up<-GOBPup[1:10,]
LungGSEAtop10down<-GOBPdown[1:10,]

uplist<-list(LungGSEAtop10up[,c(2,5,8)],ovaryGSEAtop10up[,c(2,5,8)],
             colonGSEAtop10up[,c(2,5,8)],kidneyGSEAtop10up[,c(2,5,8)])
downlist<-list(LungGSEAtop10down[,c(2,5,8)],ovaryGSEAtop10down[,c(2,5,8)],
               colonGSEAtop10down[,c(2,5,8)],kidneyGSEAtop10down[,c(2,5,8)])
upunion<-Reduce(rbind,uplist)%>%
  mutate(.,tissue=c(rep('Lung',10),rep('Ovary',10),
                    rep('Colon',10),rep('Kidney',10)))
downunion<-Reduce(rbind,downlist)%>%
  mutate(.,tissue=c(rep('Lung',10),rep('Ovary',10),
                    rep('Colon',10),rep('Kidney',10)))

pdf("GOup_intersect.pdf",height = 8,width = 8.5)
ggplot(upunion,aes(x=tissue,y=Description))+
  geom_point(aes(size=NES,color=signif(qvalue,2)))+
  theme_bw()+theme(panel.grid.minor=element_blank())+
  theme(axis.text.y = element_text(size=12,color="black"),
        axis.text.x = element_text(size=12,color="black"),
        axis.title.y = element_text(size=12,color="black"),
        axis.title.x = element_text(size=12,color="black"),)+
  scale_color_gradient(low='red',high='pink')+
  labs(color='qvalue',size='NES')+
  labs(x='',y='GOBP')+
  guides(color=guide_legend(order = 2,reverse = T),size=guide_legend(order = 1))+
  theme(text=element_text(size=12))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))
dev.off()

pdf("GOdown_intersect.pdf",height = 8,width = 8.5)
ggplot(downunion,aes(x=tissue,y=Description))+
  geom_point(aes(size=NES,color=signif(qvalue,2)))+
  theme_bw()+theme(panel.grid.minor=element_blank())+
  theme(axis.text.y = element_text(size=12,color="black"),
        axis.text.x = element_text(size=12,color="black"),
        axis.title.y = element_text(size=12,color="black"),
        axis.title.x = element_text(size=12,color="black"),)+
  scale_color_gradient(low='red',high='pink')+
  labs(color='qvalue',size='NES')+
  labs(x='',y='GOBP')+
  guides(color=guide_legend(order = 2,reverse = T),size=guide_legend(order = 1))+
  theme(text=element_text(size=12))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))
dev.off()

rm(GsGO,ErGO,genelist,GsGOres,GOBPdown,GOBPup,GSEAtop10down,GSEAtop10up,downlist,uplist,upunion,downunion)