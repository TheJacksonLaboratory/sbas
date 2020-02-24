library(limma)

library(multtest)

library(Biobase)

library(yarn)

library(edgeR)

library(ggplot2)


inc.iso.counts.mem<-read.csv('/Users/karleg/rMATs-GTEx/rmats_final.se.jc.ijc.txt',header=TRUE)


setwd('/Users/karleg/Dimorph/')

meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)


tissue.sets<-c('Heart - Left Ventricle',
               'Breast - Mammary Tissue',
               'Brain - Cortex.Brain - Frontal Cortex (BA9).Brain - Anterior cingulate cortex (BA24)',
               'Adrenal Gland'
               ,'Adipose - Subcutaneous','Muscle - Skeletal','Thyroid','Cells - Transformed fibroblasts','Artery - Aorta','Skin - Sun Exposed (Lower leg).Skin - Not Sun Exposed (Suprapubic)')

nmd.tab<-read.table('/Users/karleg/Dimorph/NMD_summary.txt',sep='\t')

an.tab<-read.table(paste0('/Users/karleg/Dimorph/fromGTF.','SE','.txt'),header=T)

an.tab<-an.tab[as.integer(nmd.tab[,2])>0,]

nmd.tab<-nmd.tab[as.integer(nmd.tab[,2])>0,]

df=data.frame(type=ifelse(as.numeric(nmd.tab[,1])>0,'NMD','NO-NMD'),length=log2(an.tab$exonEnd-an.tab$exonStart_0base))

ggplot(df, aes(x=type, y=length,color=type)) + geom_boxplot()+ylab('Log2(Exon Length)')+xlab('')+theme(axis.text=element_text(size=16))+theme(plot.background = element_rect(size=0.5,linetype="solid",color="black"))

wilcox.test((an.tab$exonEnd-an.tab$exonStart_0base)[as.numeric(nmd.tab[,1])>0],(an.tab$exonEnd-an.tab$exonStart_0base)[as.numeric(nmd.tab[,1])==0],paired=F)

mean((an.tab$exonEnd-an.tab$exonStart_0base)[as.numeric(nmd.tab[,1])>0])

mean((an.tab$exonEnd-an.tab$exonStart_0base)[as.numeric(nmd.tab[,1])==0])

df<-data.frame()

p<-list()

for (tissue.set in tissue.sets)
{
  
  if (!file.exists(paste0('/Users/karleg/Dimorph/other/',tissue.set,'se.txt')))
    
    next
  
  
  meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)
  
  print(tissue.set)
  
  meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)
  
  inc.iso.counts<-inc.iso.counts.mem
  
  meta.data$Sample_Name_s<-gsub('-','\\.',meta.data$Sample_Name_s)
  
  meta.data<-meta.data[meta.data$Sample_Name_s %in% colnames(inc.iso.counts),]
  
  rownames(inc.iso.counts)<-inc.iso.counts$ID
  
  inc.iso.counts<-inc.iso.counts[,-1]
  
  inc.iso.counts<-inc.iso.counts[,order(match(colnames(inc.iso.counts),meta.data$Sample_Name_s))]
  
  inc.iso.counts<-inc.iso.counts[,!grepl('11ILO',meta.data$Sample_Name_s) & (meta.data$body_site_s %in% unlist(strsplit(tissue.set,split='\\.')))]
  
  ds.tab<-read.table(paste0('/Users/karleg/Dimorph/other/',paste(tissue.set,collapse = '.'),'se.txt'))
  
  ds.tab<-ds.tab[abs(ds.tab$logFC)>=log2(1.5) & ds.tab$adj.P.Val<=0.05,]
  
  mean(rowMeans(inc.iso.counts)[rownames(inc.iso.counts) %in% rownames(ds.tab)])
  
  ds.rowmeans<-rowMeans(inc.iso.counts)[rownames(inc.iso.counts) %in% rownames(ds.tab)]
  
  nmd.tab<-read.table('/Users/karleg/Dimorph/NMD_summary.txt',sep='\t')
  
  hbm.tab<-read.table('/Users/karleg/Dimorph/summary_hbm.txt',header=T,sep='\t')
  
  hbm.tab<-hbm.tab[hbm.tab$Tissue==tissue.set & hbm.tab$Dimorph=='Yes' & hbm.tab$Sig..Gene.Expression!=0,]
  
  hbm.tab<-hbm.tab[which(nmd.tab$num.transcripts[hbm.tab$Event]>0),]
  
  df<-rbind(df,data.frame(type=c(ifelse(as.numeric(nmd.tab[hbm.tab$Event,1])>0,'NMD','NO-NMD')),log2counts=c(log2(rowMeans(inc.iso.counts[hbm.tab$Event,])+0.5))))
    
}
  
  

pval<-wilcox.test(df$log2counts[df$type=='NMD'],df$log2counts[df$type=='NO-NMD'],paired=F)$p.value
  
p <- ggplot(df, aes(x=type, y=log2counts,color=type)) + geom_boxplot()+theme(plot.title = element_text(hjust = 0.5))+ylab('log2(counts+0.5)')+xlab('')+
    theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"))+ theme(legend.position="none")+
    annotate("text", x=1.5, y=10, label= paste0("p-value: ",round(pval,-log10(pval)+3))) 
p  

save.image('/Users/karleg/Dimorph/RDATA/figure4d.RData')


#Before running the following, use the Session menu to set working directory to source file location
load('figure4d.RData')
library(ggplot2)
library(ggsci)


pn<-ggplot(df, aes(x=type, y=log2counts,fill=type))  + geom_boxplot(notch=TRUE) + theme_minimal() + scale_fill_npg() 
pn<- pn + theme(axis.text = element_text(size=24, hjust=0.5),
          axis.title.x=element_blank(),
          axis.title.y = element_text(face="bold",  size=22),
          axis.text.y = element_text(size=36),
          legend.position="none")   +ylab(expression(log[2]~(counts+0.5)))
my.pval = "p==1.55 %*% 10 ^ -10"  ##
pn <- pn + annotate("text", x=1.5, y=10, label= as.character(my.pval),parse=TRUE, fontface="bold", size=10) 
pn
