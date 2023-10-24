library(ggvenn)
#load('FFPE_object.RData')
###paired DEG analysis for Lung tissue###
countdata<-Lung_code_count
dds<- DESeqDataSetFromMatrix(countData = as.matrix(countdata),
                             colData = Lung_paired_meta,
                             design = ~sampleid+FF)
relevel(dds$FF, ref = "FF")->dds$condition 
dds <- DESeq(dds)
Lung_deg<-data.frame(results(dds))
Lung_deg %>%
  mutate(group = case_when(
    baseMean > 1 & log2FoldChange >= 1 & padj < 0.05 ~ "up",
    baseMean > 1 & log2FoldChange <= -1 & padj < 0.05 ~ "down",
    TRUE ~ "not change"
  ),Geneid = rownames(Lung_deg))->Lung_deg

###paired DEG analysis for BPV program data###
target<-BPVpaired_meta$Body.Site=='COLON'
condition = factor(BPVpaired_meta$FF[target])
sampletype = factor(BPVpaired_meta$Body.Site[target])
sampleid <- factor(BPVpaired_meta$sampleid[target])
coldata<-data.frame(row.names = BPVpaired_meta$Run[target],condition,sampleid,sampletype)
countdata<-BPV_count%>%
  select(.,BPVpaired_meta$Run[target])
dds<- DESeqDataSetFromMatrix(countData = as.matrix(countdata),
                             colData = coldata,
                             design = ~sampleid+condition)
relevel(dds$condition, ref = "FF")->dds$condition 
dds <- DESeq(dds)
colon_deg<-data.frame(results(dds))
colon_deg %>%
  mutate(group = case_when(
    baseMean > 1 & log2FoldChange >= 1 & padj < 0.05 ~ "up",
    baseMean > 1 & log2FoldChange <= -1 & padj < 0.05 ~ "down",
    TRUE ~ "not change"
  ),Geneid = rownames(colon_deg))->colon_deg


###deg venn plot####
filter(colon_deg,padj<0.05&baseMean>1&log2FoldChange>=1)%>%
  rownames()->tmp
filter(Lung_deg,padj<0.05&baseMean>1&log2FoldChange>=1)%>%
  rownames()->tmp1
filter(kidney_deg,padj<0.05&baseMean>1&log2FoldChange>=1)%>%
  rownames()->tmp2
filter(ovary_deg,padj<0.05&baseMean>1&log2FoldChange>=1)%>%
  rownames()->tmp3
pdf("FFPEup_venn.pdf",width = 6,height = 5)
ggvenn(list('Colon'=tmp,
            'Lung'=tmp1,
            'Kidney'=tmp2,
            'Ovary'=tmp3
            ),
       stroke_color = "white",
       fill_color = c("blue", "orange", "green", "red"),
       set_name_color = c("blue", "orange", "green", "red"))+
  theme(plot.title = element_text(vjust = 0,hjust = 0.5,size = 20))
dev.off()

filter(colon_deg,padj<0.05&baseMean>1&log2FoldChange<=-1)%>%
  rownames()->tmp
filter(Lung_deg,padj<0.05&baseMean>1&log2FoldChange<=-1)%>%
  rownames()->tmp1
filter(kidney_deg,padj<0.05&baseMean>1&log2FoldChange<=-1)%>%
  rownames()->tmp2
filter(ovary_deg,padj<0.05&baseMean>1&log2FoldChange<=-1)%>%
  rownames()->tmp3
pdf("FFPEdn_venn.pdf",width = 6,height = 5)
ggvenn(list('Colon'=tmp,
            'Lung'=tmp1,
            'Kidney'=tmp2,
            'Ovary'=tmp3
),
stroke_color = "white",
fill_color = c("blue", "orange", "green", "red"),
set_name_color = c("blue", "orange", "green", "red"))+
  #labs(title="downregulated genes in FFPE",)+
  theme(plot.title = element_text(vjust = 0,hjust = 0.5,size = 20))
dev.off()