library(cowplot)
library(patchwork)
#load('FFPE_object.RData')
####import dataset#############
path <-"/rseqc_output/summary"
filenames<-dir(path)
filenames<-filenames[grep('^res',filenames)]
filepath <-sapply(filenames,function(x){
  paste(path,x,sep = '/')})
read_distribution_file1<-lapply(filepath, function(x){
  read.csv(x, sep="")})
path <-"/rseqc_output/summary"
filenames<-dir(path)
filepath <-sapply(filenames,function(x){
  paste(path,x,sep = '/')})
read_distribution_file2<-lapply(filepath, function(x){
  read.csv(x, sep="")})

for (i in 2:length(read_distribution_file1)) {
  if(i==2){
    read_distribution1<-merge(read_distribution_file1[[1]][1:4,c(1,3)],read_distribution_file1[[2]][1:4,c(1,3)],by='Group')}
  else
    read_distribution1<-merge(read_distribution1,read_distribution_file1[[i]][1:4,c(1,3)],by='Group')
}
colnames(read_distribution1)[-1]<-c('L1_FFPE','L4_FF','L5_FF','L6_FF','L2_FFPE','L3_FFPE',
                                    'L4_FFPE','L5_FFPE','L6_FFPE','L1_FF','L2_FF','L3_FF')
read_distribution1[,-1]<-apply(read_distribution1[,-1],2,
                                function(x){x/sum(x)})
read_distribution1%>%gather(.,key = "sample",value = "proportion",-Group)%>%
  data.frame()->read_distribution1_narrow
read_distribution1_narrow$sampletype<-sapply(read_distribution1_narrow$sample,function(x){
  unlist(strsplit(x,split = "_"))[2]})
read_distribution1_narrow$sampleid<-sapply(read_distribution1_narrow$sample,function(x){
  unlist(strsplit(x,split = "_"))[1]})

for (i in 2:length(read_distribution_file2)) {
  if(i==2){
    read_distribution2<-merge(read_distribution_file2[[1]][1:4,c(1,3)],read_distribution_file2[[2]][1:4,c(1,3)],by='Group')}
  else
    read_distribution2<-merge(read_distribution2,read_distribution_file2[[i]][1:4,c(1,3)],by='Group')
}
colnames(read_distribution2)[-1]<-c('L7_FF','L7_40.5','L7_406','L7_448','L7_c0.5','L7_c06',
                                    'L7_c48','L3.1_FF','L3.1_FFPE','L4.1_FF','L4.1_FFPE','L7_section')
read_distribution2[,-1]<-apply(read_distribution2[,-1],2,
                               function(x){x/sum(x)})
read_distribution2%>%gather(.,key = "sample",value = "proportion",-Group)%>%
  data.frame()->read_distribution2_narrow
read_distribution2_narrow$sampletype<-sapply(read_distribution2_narrow$sample,function(x){
  unlist(strsplit(x,split = "_"))[2]})
read_distribution2_narrow$sampleid<-sapply(read_distribution2_narrow$sample,function(x){
  unlist(strsplit(x,split = "_"))[1]})

read_distribution1_narrow$Group<-gsub('_',' ',read_distribution1_narrow$Group)
read_distribution2_narrow$Group<-gsub('_',' ',read_distribution2_narrow$Group)
a<-ggplot(rbind(read_distribution1_narrow,read_distribution2_narrow[29:44,]), aes( x = sampletype,y=proportion,fill = Group))+
  geom_col(position = 'stack', width = 0.8)+
  facet_wrap(~sampleid,nrow = 1,strip.position = 'bottom')+
  theme(axis.text.x=element_text(colour = "black",size=14,angle = 45,hjust = 1),
        axis.text.y=element_text(colour = "black",size=14),
        axis.title.y=element_text(colour = "black",size=14),
        legend.text = element_text(colour = "black",size=14),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal")+
  xlab('')+ylab('Proportion')

tmp<-read_distribution2_narrow[-29:-44,]
tmp$sampletype<-factor(tmp$sampletype,levels = unique(tmp$sampletype))
b<-ggplot(tmp, aes( x = sampletype,y=proportion,fill = Group))+
  geom_col(position = 'stack', width = 0.8)+
  facet_wrap(~sampleid,nrow = 1,strip.position = 'bottom')+
  theme(axis.text.x=element_text(colour = "black",size=14,angle = 45,hjust = 1),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(colour = "black",size=14),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal")+
  xlab('')+ylab('Proportion')

pdf('read_distribution.pdf',height = 4.7,width = 9.5)
wrap_plots(A=a,B=b,nrow = 1,design = c('AAB'))+plot_layout(guides = 'collect')&
  theme(legend.position='bottom')
dev.off()

