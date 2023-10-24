#to assess the ultility of FFPE
library(grDevices)
library(RColorBrewer)
library(FactoMineR)#PCA function
library(factoextra)#visualization
library(tidyr)
library(dplyr)
library(DESeq2)
library(ggrepel)
library(pheatmap)
library(eulerr)
#load('FFPE_object.RData')
####filter genes before PCA####
df<-vst(as.matrix(FFPE1count))
df<-df[apply(df,1,mad)>0.2,]
gene.pca<-PCA(t(df), ncp = 5, scale.unit = T, graph = FALSE)
####hierarchical cluster visualization
sampleTree = hclust(dist(gene.pca$ind$coord))
pdf(file = "hclust.pdf",width = 4,height = 4)
fviz_dend(sampleTree,k=4,rect =T,rect_fill = T,main = 'Cluster dendrogram')
dev.off()

d<-dist(gene.pca$ind$coord)
pheatmap(as.matrix(d),treeheight_row = 0,
         display_numbers = F,fontsize = 15,
         filename = "hclust1.pdf")


Grp<-c(rep(c("Tumor","Peritumor"),6))
Lbl<-c("FF.1T","FF.1P","FF.2T","FF.2P","FF.3T","FF.3P",
       "FFPE.1T","FFPE.1P","FFPE.2T","FFPE.2P","FFPE.3T","FFPE.3P")
##visualization
#`protein_coding_104` content ensembl_id/hgnc_symbol information of 19349 gene
gene.pca <- PCA(t(df[,1:12]), ncp = 5, scale.unit = T, graph = FALSE)
plot(gene.pca)
pca_sample <- data.frame(gene.pca$ind$coord[,c(1,2)])
pca_sample$group<-as.factor(Grp[1:12])
pca_sample$label<-as.factor(Lbl[1:12])
pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )
pdf(file = "PCA_FF_FFPE.pdf",width = 4.5,height = 4.5)
ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = group), size = 2) + 
  geom_text_repel(aes(label=label),size=4,max.overlaps = 20)+
  theme(legend.text = element_text(size=14),axis.text.x = element_text(size=12,color = 'black'),
        axis.text.y = element_text(size=12,color = 'black'),axis.title.x =  element_text(size=14)
        ,axis.title.y =  element_text(size=14),panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),legend.position='bottom') +  #去除背景和网格线
  labs(title="",x =  paste('PC1:', pca_eig1, '%'), y = paste('PC2:', pca_eig2, '%'), color = '')
dev.off()

rm(Grp,Lbl,gene.pca,pca_sample,PC1,PCAplot,pca_eig1,pca_eig2)

##dimension reduction
gene.pca <- PCA(t(df[,7:12]), ncp = 5, scale.unit = T, graph = FALSE)
fviz_pca_ind(gene.pca,col.ind = "contrib")#fviz_pca_ind(gene.pca)
Lung_pca_contribution_FFPE<-gene.pca$svd$V
gene.pca <- PCA(t(df[,1:6]), ncp = 5, scale.unit = T, graph = FALSE)
fviz_pca_ind(gene.pca,col.ind = "contrib")#fviz_pca_ind(gene.pca)
Lung_pca_contribution<-gene.pca$svd$V

PC_top_cor<-data.frame(abs(cor(Lung_pca_contribution,Lung_pca_contribution_FFPE)))
colnames(PC_top_cor)<-c(paste0("Dim.",1:5))
rownames(PC_top_cor)<-c(paste0("Dim.",1:5))

####visualize the correlation of top 5 PC
tmp<-PC_top_cor%>%
  rownames_to_column(.,var = "FF")%>%
  gather(.,value = "Pearson r",key = "FFPE",-FF)
tmp$FF<-str_replace_all(tmp$FF,'Dim','PC')
tmp$FFPE<-str_replace_all(tmp$FFPE,'Dim','PC')
tmp$FFPE<-factor(tmp$FFPE,levels = unique(tmp$FFPE)[5:1])

pdf(file = "top5PC.pdf",width = 6,height = 5)
ggplot(tmp, aes(x = FF, y = FFPE)) +
  geom_tile(aes(fill = `Pearson r`)) +
  geom_text(aes(label = round(`Pearson r`, 2),size=2))+
  scale_fill_gradient(low = "white", high = "red") +
  guides(size="none")+
  theme(legend.text = element_text(size=12,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size=12,color='black'),
        axis.title.x =  element_text(size=14,color='black')
        ,axis.title.y =  element_text(size=14),panel.grid = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.key = element_rect(fill = 'transparent'))
dev.off()

####deseq2_deg#####
countdata<-as.matrix(FFPE1count)
dds<- DESeqDataSetFromMatrix(countData = countdata[,c(1:6)],
                             colData = FFPE1metadata[c(1:6),],
                             design = ~pair+condition)
relevel(dds$condition, ref = "P")->dds$condition 
dds <- DESeq(dds)
resultsNames(dds)
summary(results(dds))
P_A_deg_FF_pair<-data.frame(results(dds,name = "condition_A_vs_P" ))
P_A_deg_FF_pair %>%
  mutate(group = case_when(
    baseMean > 10 & log2FoldChange >= 1 & padj < 0.05 ~ "up",
    baseMean > 10 & log2FoldChange <= -1 & padj < 0.05 ~ "down",
    TRUE ~ "not significant"
  ),Geneid = protein_coding$Geneid)->P_A_deg_FF_pair

dds<- DESeqDataSetFromMatrix(countData = countdata[,c(7:12)],
                             colData = FFPE1metadata[c(7:12),],
                             design = ~pair+condition)
relevel(dds$condition, ref = "P")->dds$condition 
dds <- DESeq(dds)
resultsNames(dds)
summary(results(dds))
P_A_deg_FFPE_pair<-data.frame(results(dds,name = "condition_A_vs_P" ))
P_A_deg_FFPE_pair %>%
  mutate(group = case_when(
    baseMean > 10 & log2FoldChange >= 1 & padj < 0.05 ~ "up",
    baseMean > 10 & log2FoldChange <= -1 & padj < 0.05 ~ "down",
    TRUE ~ "not significant"
  ),Geneid = protein_coding$Geneid)->P_A_deg_FFPE_pair

####DEG result consistency
library(ggvenn)
####deg venn plot####
P_A_deg_FF_pair[,]%>%
filter(.,padj<0.05&baseMean>1&log2FoldChange> 1)%>%
  rownames()->tmp
P_A_deg_FFPE_pair[,]%>%
filter(.,padj<0.05&baseMean>1&log2FoldChange> 1)%>%
  rownames()->tmp1

dat1<-c(" " = length(setdiff(tmp,tmp1)), 
       "  " = length(setdiff(tmp1,tmp)), 
       " &  " = length(intersect(tmp1,tmp)))
pdf(file = "up_overlap.pdf",width = 6.5,height = 4.5)
plot(euler(dat1),
     fills = list(fill=c("#E41A1C","#1E90FF",'#7B68EE'),
                  alpha=0.5),
     quantities = list(c(length(setdiff(tmp,tmp1)),length(intersect(tmp1,tmp)),length(setdiff(tmp1,tmp))),
                       col = 'black',
                       cex =4),
     lebels = list(col = 'blank'))
dev.off()

####kegg富集分析一致性####
upunion<-rbind(kegg_up_FF_top300@result[1:5,c(2,3,7,9)],kegg_up_FFPE_top300@result[1:5,c(2,3,7,9)])%>%
  mutate(group=c(rep('FF',5),rep('FFPE',5)),
         GeneRatio=signif(c(Count[1:5]/123,Count[6:10]/135),digits = 2))
downunion<-rbind(kegg_dn_FF_top300@result[1:5,c(2,3,7,9)],kegg_dn_FFPE_top300@result[1:5,c(2,3,7,9)])%>%
  mutate(group=c(rep('FF',5),rep('FFPE',5)),
         GeneRatio=signif(c(Count[1:5]/152,Count[6:10]/154),digits = 2))

pdf("/home/liny/R/labworkbase/FFPE_final/FFPE_final_plot/plot/utility/upkegg_intersect_top300.pdf",height = 5.5,width = 6)
ggplot(upunion,aes(x=group,y=Description))+
  geom_point(aes(size=GeneRatio,color=signif(qvalue,2)))+
  theme_bw()+theme(panel.grid.minor=element_blank())+
  theme(axis.text.y = element_text(size=14,color="black"),
        axis.text.x = element_text(size=14,color="black"),
        axis.title.y = element_text(size=14,color="black"),
        axis.title.x = element_text(size=12,color="black"),
        axis.title = element_text(size=14,color = 'black'))+
  scale_color_gradient(low='red',high='pink')+
  labs(color='qvalue',size='GeneRatio')+
  labs(x='',y='KEGG',title = "upregulated")+
  guides(color=guide_legend(order = 2,reverse = T),size=guide_legend(order = 1))+
  theme(text=element_text(size=14))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=30))
dev.off()

pdf("/home/liny/R/labworkbase/FFPE_final/FFPE_final_plot/plot/utility/dnkegg_intersect_top300.pdf",height = 5.5,width = 6)
ggplot(downunion,aes(x=group,y=Description))+
  geom_point(aes(size=GeneRatio,color=signif(qvalue,2)))+
  theme_bw()+theme(panel.grid.minor=element_blank())+
  theme(axis.text.y = element_text(size=14,color="black"),
        axis.text.x = element_text(size=14,color="black"),
        axis.title.y = element_text(size=14,color="black"),
        axis.title.x = element_text(size=14,color="black"),
        axis.title = element_text(size=14,color = 'black'))+
  scale_color_gradient(low='red',high='pink')+
  labs(color='qvalue',size='GeneRatio')+
  labs(x='',y='KEGG',title = "downregulated")+
  guides(color=guide_legend(order = 2,reverse = T),size=guide_legend(order = 1))+
  theme(text=element_text(size=14))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=30))
dev.off()

####consistency of immune and pan-cancer gene set enrichment analyses####
immune_list1 <- read.delim2("immuport_GeneList.txt")
immune_list_name<-unique(immune_list1$Category)
immune_list <-list()
for (item in immune_list_name) {
  immune_list[[item]]<-immune_list1$Symbol[immune_list1$Category==item]
}
for (i in 1:17) {
  if(i==1){immune_geneset<-list()}
  immune_geneset[[i]]<-GSEABase::GeneSet(unique(immune_list[[i]]),setName=immune_list_name[i],geneIdType=SymbolIdentifier())
}
immune_geneset <- GSEABase::GeneSetCollection(immune_geneset)

#cancer gene list
nanostring_cancerfunction <- read.csv("nanostring_cancerfunction.csv")
nanostring_immune <- read.csv("nanostring_immune.csv")
nanostring_cancerfunction_narrow<-tidyr::gather(nanostring_cancerfunction,key='pathway',value='ifin',-Gene)%>%filter(ifin=='+')
nanostring_cancerfunction_narrow$pathway<-str_replace_all(nanostring_cancerfunction_narrow$pathway,pattern = '\\.',replacement = ' ')
nanostring_immune_narrow<-tidyr::gather(nanostring_immune,key='pathway',value='ifin',-Gene)%>%filter(ifin=='+')
nanostring_immune_narrow$pathway<-str_replace_all(nanostring_immune_narrow$pathway,pattern = '\\.',replacement = ' ')
#immune gene list
cancerimmune_list_name<-unique(nanostring_immune_narrow$pathway)
cancerimmune_list <-list()
for (item in cancerimmune_list_name) {
  cancerimmune_list[[item]]<-nanostring_immune_narrow$Gene[nanostring_immune_narrow$pathway==item]
}
for (i in 1:length(cancerimmune_list_name)) {
  if(i==1){cancerimmune_geneset<-list()}
  cancerimmune_geneset[[i]]<-GSEABase::GeneSet(unique(cancerimmune_list[[i]]),setName=cancerimmune_list_name[i],geneIdType=SymbolIdentifier())
}
cancerimmune_geneset <- GSEABase::GeneSetCollection(cancerimmune_geneset)
#pan-cancer pathway
cancer_list_name<-unique(nanostring_cancerfunction_narrow$pathway)
cancer_list <-list()
for (item in cancer_list_name) {
  cancer_list[[item]]<-nanostring_cancerfunction_narrow$Gene[nanostring_cancerfunction_narrow$pathway==item]
}
for (i in 1:length(cancer_list_name)) {
  if(i==1){cancer_geneset<-list()}
  cancer_geneset[[i]]<-GSEABase::GeneSet(unique(cancer_list[[i]]),setName=cancer_list_name[i],geneIdType=SymbolIdentifier())
}
cancer_geneset <- GSEABase::GeneSetCollection(cancer_geneset)

###ORA分析（不太合理）
P_A_deg_FF_pair<-P_A_deg_FF_pair%>%
  mutate("Pancancer_panel"=rownames(.)%in%pancancergenelist$Gene.Name,
         "Antigen_Processing_and_Presentation" = rownames(.)%in%immune_list[["Antigen_Processing_and_Presentation"]],
         "Antimicrobials" = rownames(.)%in%immune_list[["Antimicrobials"]],
         "BCRSignalingPathway" = rownames(.)%in%immune_list[["BCRSignalingPathway"]],
         "Chemokines" = rownames(.)%in%immune_list[["Chemokines"]],
         "Cytokines" = rownames(.)%in%immune_list[["Cytokines"]],
         "Cytokine_Receptors" = rownames(.)%in%immune_list[["Cytokine_Receptors"]],
         "Chemokine_Receptors" = rownames(.)%in%immune_list[["Chemokine_Receptors"]],
         "Interferons" = rownames(.)%in%immune_list[["Interferons"]],
         "Interferon_Receptor" = rownames(.)%in%immune_list[["Interferon_Receptor"]],
         "Interleukins" = rownames(.)%in%immune_list[["Interleukins"]],
         "Interleukins_Receptor" = rownames(.)%in%immune_list[["Interleukins_Receptor"]],
         "NaturalKiller_Cell_Cytotoxicity" = rownames(.)%in%immune_list[["NaturalKiller_Cell_Cytotoxicity"]],
         "TCRsignalingPathway" = rownames(.)%in%immune_list[["TCRsignalingPathway"]],
         "TGFb_Family_Member" = rownames(.)%in%immune_list[["TGFb_Family_Member"]],
         "TGFb_Family_Member_Receptor" = rownames(.)%in%immune_list[["TGFb_Family_Member_Receptor"]],
         "TNF_Family_Members" = rownames(.)%in%immune_list[["TNF_Family_Members"]],
         "TNF_Family_Members_Receptors" = rownames(.)%in%immune_list[["TNF_Family_Members_Receptors"]],
         "significant" = group!="not significant"
         )
P_A_deg_FFPE_pair<-P_A_deg_FFPE_pair%>%
  mutate("Pancancer_panel"=rownames(.)%in%pancancergenelist$Gene.Name,
         "Antigen_Processing_and_Presentation" = rownames(.)%in%immune_list[["Antigen_Processing_and_Presentation"]],
         "Antimicrobials" = rownames(.)%in%immune_list[["Antimicrobials"]],
         "BCRSignalingPathway" = rownames(.)%in%immune_list[["BCRSignalingPathway"]],
         "Chemokines" = rownames(.)%in%immune_list[["Chemokines"]],
         "Cytokines" = rownames(.)%in%immune_list[["Cytokines"]],
         "Cytokine_Receptors" = rownames(.)%in%immune_list[["Cytokine_Receptors"]],
         "Chemokine_Receptors" = rownames(.)%in%immune_list[["Chemokine_Receptors"]],
         "Interferons" = rownames(.)%in%immune_list[["Interferons"]],
         "Interferon_Receptor" = rownames(.)%in%immune_list[["Interferon_Receptor"]],
         "Interleukins" = rownames(.)%in%immune_list[["Interleukins"]],
         "Interleukins_Receptor" = rownames(.)%in%immune_list[["Interleukins_Receptor"]],
         "NaturalKiller_Cell_Cytotoxicity" = rownames(.)%in%immune_list[["NaturalKiller_Cell_Cytotoxicity"]],
         "TCRsignalingPathway" = rownames(.)%in%immune_list[["TCRsignalingPathway"]],
         "TGFb_Family_Member" = rownames(.)%in%immune_list[["TGFb_Family_Member"]],
         "TGFb_Family_Member_Receptor" = rownames(.)%in%immune_list[["TGFb_Family_Member_Receptor"]],
         "TNF_Family_Members" = rownames(.)%in%immune_list[["TNF_Family_Members"]],
         "TNF_Family_Members_Receptors" = rownames(.)%in%immune_list[["TNF_Family_Members_Receptors"]],
         "significant" = group!="not significant"
  )

#GSEA
genelist<-P_A_deg_FFPE_pair[P_A_deg_FFPE_pair$baseMean>10,2]
names(genelist)<-rownames(P_A_deg_FFPE_pair)[P_A_deg_FFPE_pair$baseMean>10]
genelist<-sort(genelist,decreasing = T)

tmp<-GSEA(geneList = genelist,pvalueCutoff = 0.05,TERM2GENE = rbind(nanostring_cancerfunction_narrow[,c(2,1)],nanostring_immune_narrow[,c(2,1)]))

genelist<-P_A_deg_FF_pair[P_A_deg_FF_pair$baseMean>10,2]
names(genelist)<-rownames(P_A_deg_FF_pair)[P_A_deg_FF_pair$baseMean>10]
genelist<-sort(genelist,decreasing = T)

tmp1<-GSEA(geneList = genelist,pvalueCutoff = 0.05,TERM2GENE = rbind(nanostring_cancerfunction_narrow[,c(2,1)],nanostring_immune_narrow[,c(2,1)]))

gse_union<-rbind(tmp1@result,tmp@result[tmp@result$NES>0,])
gse_union$Description<-as.factor(gse_union$Description)
gse_union$group<-factor(c(rep("FF",8),rep("FFPE",10)),levels=c("FF","FFPE"))

pdf("GSEA_cancer_overlap.pdf",height = 5.5,width = 4.5)
ggplot(gse_union,aes(x=group,y=Description))+
  geom_point(aes(size=NES,color=signif(p.adjust,2)))+
  theme_bw()+theme(panel.grid.minor=element_blank())+
  theme(axis.text.y = element_text(size=12,color="black"),
        axis.text.x = element_text(size=12,color="black"),
        axis.title.y = element_text(size=12,color="black"),
        axis.title.x = element_text(size=12,color="black"),
        legend.text = element_text(size=12,color="black"),
        legend.title = element_text(size=12,color="black"))+
  scale_color_gradient(low='red',high='pink')+
  labs(color='p.adjust',size='NES')+
  labs(x='',y='Annotations')+
  guides(color=guide_legend(order = 2,reverse = T),size=guide_legend(order = 1))+
  theme(text=element_text(size=12))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=25))
dev.off()
