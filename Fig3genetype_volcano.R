library(ggplot2)
library(ggpubr)
library(ggrepel)
#load('FFPE_object.RData')
taste_pattern='^TAS[12]R*'
Olfactory_pattern = '^OR[1-9][A-Z]*'

to_check<-Reduce(f = union,x = list(grep(colon_deg$hgnc_symbol,pattern = taste_pattern),
                          grep(colon_deg$hgnc_symbol,pattern = Olfactory_pattern)))
Lung_deg$genetype<-'Other'
Lung_deg$genetype[grep(colon_deg$hgnc_symbol,pattern = taste_pattern)]<-'Taste'
Lung_deg$genetype[grep(colon_deg$hgnc_symbol,pattern = Olfactory_pattern)]<-'Olfactory'

data<-Lung_deg[to_check,]
data$genetype<-factor(data$genetype,levels=c('Taste','Olfactory'))
pdf("Lung_volcano.pdf",width = 7,height = 5)
ggplot(data[data$baseMean>1,], aes(x = log2FoldChange, y = -log10(padj),color=genetype)) +
  geom_point(data=Lung_deg[Lung_deg$genetype=='Other',],mapping = aes(x = log2FoldChange, y = -log10(padj)),color='grey',alpha=0.2,size=2)+
  geom_point(alpha=0.8,size=2) +
  labs(title = "Differential expressed genes in lung", x = "log2FoldChange", y = "-log[10](padj)")+
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12,color = 'black'),
        axis.title.x = element_text(size=12,color = 'black'),
        axis.text.y = element_text(size=12,color = 'black'),
        axis.title.y = element_text(size=12,color = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(size=12,color = 'black'),
        legend.position = 'bottom')
dev.off()

pdf("Lung_violin.pdf",width = 4,height = 5)
ggplot(data[data$baseMean>1,],aes(y=log2FoldChange,x=genetype,fill=genetype))+
  geom_violin(trim = FALSE, scale = "width", alpha = 1)+
  geom_boxplot(width=0.3,fill='white')+
  geom_hline(yintercept = c(-1,1), linetype = "dashed")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12,color = 'black'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12,color = 'black'),
        axis.title.y = element_text(size=12,color = 'black'),
        legend.position = 'none')
dev.off()

colon_deg$genetype<-'Other'
colon_deg$genetype[grep(colon_deg$hgnc_symbol,pattern = taste_pattern)]<-'Taste'
colon_deg$genetype[grep(colon_deg$hgnc_symbol,pattern = Olfactory_pattern)]<-'Olfactory'

data<-colon_deg[to_check,]
data$genetype<-factor(data$genetype,levels=c('Taste','Olfactory'))
pdf("colon_volcano.pdf",width = 7,height = 5)
ggplot(data[data$baseMean>1,], aes(x = log2FoldChange, y = -log10(padj),color=genetype)) +
  geom_point(data=colon_deg[colon_deg$genetype=='Other',],mapping = aes(x = log2FoldChange, y = -log10(padj)),color='grey',alpha=0.2,size=2)+
  geom_point(alpha=0.8,size=2) +
  labs(title = "Differential expressed genes in colon", x = "log2FoldChange", y = "-log[10](padj)")+
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12,color = 'black'),
        axis.title.x = element_text(size=12,color = 'black'),
        axis.text.y = element_text(size=12,color = 'black'),
        axis.title.y = element_text(size=12,color = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(size=12,color = 'black'),
        legend.position = 'bottom')+
  xlim(-5,5)
dev.off()

pdf("colon_violin.pdf",width = 4,height = 5)
ggplot(data[data$baseMean>1,],aes(y=log2FoldChange,x=genetype,fill=genetype))+
  geom_violin(trim = FALSE, scale = "width", alpha = 1)+
  geom_boxplot(width=0.3,fill='white')+
  geom_hline(yintercept = c(-1,1), linetype = "dashed")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12,color = 'black'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12,color = 'black'),
        axis.title.y = element_text(size=12,color = 'black'),
        legend.position = 'none')
dev.off()

kidney_deg$genetype<-'Other'
kidney_deg$genetype[grep(colon_deg$hgnc_symbol,pattern = taste_pattern)]<-'Taste'
kidney_deg$genetype[grep(colon_deg$hgnc_symbol,pattern = Olfactory_pattern)]<-'Olfactory'

data<-kidney_deg[to_check,]
data$genetype<-factor(data$genetype,levels=c('Taste','Olfactory'))
pdf("kidney_volcano.pdf",width = 7,height = 5)
ggplot(data[data$baseMean>1,], aes(x = log2FoldChange, y = -log10(padj),color=genetype)) +
  geom_point(data=kidney_deg[kidney_deg$genetype=='Other',],mapping = aes(x = log2FoldChange, y = -log10(padj)),color='grey',alpha=0.2,size=2)+
  geom_point(alpha=0.8,size=2) +
  labs(title = "Differential expressed genes in kidney", x = "log2FoldChange", y = "-log[10](padj)")+
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12,color = 'black'),
        axis.title.x = element_text(size=12,color = 'black'),
        axis.text.y = element_text(size=12,color = 'black'),
        axis.title.y = element_text(size=12,color = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(size=12,color = 'black'),
        legend.position = 'bottom')+
  xlim(-5,5)
dev.off()

pdf("kidney_violin.pdf",width = 4,height = 5)
ggplot(data[data$baseMean>1,],aes(y=log2FoldChange,x=genetype,fill=genetype))+
  geom_violin(trim = FALSE, scale = "width", alpha = 1)+
  geom_boxplot(width=0.3,fill='white')+
  geom_hline(yintercept = c(-1,1), linetype = "dashed")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12,color = 'black'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12,color = 'black'),
        axis.title.y = element_text(size=12,color = 'black'),
        legend.position = 'none')
dev.off()

ovary_deg$genetype<-'Other'
ovary_deg$genetype[grep(colon_deg$hgnc_symbol,pattern = taste_pattern)]<-'Taste'
ovary_deg$genetype[grep(colon_deg$hgnc_symbol,pattern = Olfactory_pattern)]<-'Olfactory'

data<-ovary_deg[to_check,]
data$genetype<-factor(data$genetype,levels=c('Taste','Olfactory'))
pdf("ovary_volcano.pdf",width = 7,height = 5)
ggplot(data[data$baseMean>1,], aes(x = log2FoldChange, y = -log10(padj),color=genetype)) +
  geom_point(data=ovary_deg[ovary_deg$genetype=='Other',],mapping = aes(x = log2FoldChange, y = -log10(padj)),color='grey',alpha=0.2,size=2)+
  geom_point(alpha=0.8,size=2) +
  labs(title = "Differential expressed genes in ovary", x = "log2FoldChange", y = "-log[10](padj)")+
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12,color = 'black'),
        axis.title.x = element_text(size=12,color = 'black'),
        axis.text.y = element_text(size=12,color = 'black'),
        axis.title.y = element_text(size=12,color = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(size=12,color = 'black'),
        legend.position = 'bottom')+
  xlim(-5,5)
dev.off()

pdf("ovary_violin.pdf",width = 4,height = 5)
ggplot(data[data$baseMean>1,],aes(y=log2FoldChange,x=genetype,fill=genetype))+
  geom_violin(trim = FALSE, scale = "width", alpha = 1)+
  geom_boxplot(width=0.3,fill='white')+
  geom_hline(yintercept = c(-1,1), linetype = "dashed")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12,color = 'black'),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12,color = 'black'),
        axis.title.y = element_text(size=12,color = 'black'),
        legend.position = 'none')
dev.off()

