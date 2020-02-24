

meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

tissues<-unique(meta.data$body_site_s)

all.genes<-read.table('/Users/karleg/Dimorph/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct',sep='\t',header=T,skip=2,colClasses = c(rep("character", 2), rep("NULL", 11688)))

all.genes<-all.genes[!duplicated(all.genes[,'Description']),'Description']

fc.exp.mat<-matrix(rep(0), length(all.genes)*length(tissues),nrow = length(all.genes),ncol=length(tissues))

colnames(fc.exp.mat)<-tissues

remove=c()

for (tissue in tissues)
{
  
  if (!file.exists(paste0('/Users/karleg/Dimorph/gene_expression/all_genes/DE_result_',tissue,'.txt')))
  {
    remove<-c(remove,which(colnames(fc.exp.mat)==tissue))
    
  }else{
    
    tab<-read.table(paste0('/Users/karleg/Dimorph/gene_expression/all_genes/DE_result_',tissue,'.txt'))
    
    tab<-tab[order(match(rownames(tab),all.genes)),]
    
    fc.exp.mat[which(all.genes %in% rownames(tab)),which(colnames(fc.exp.mat)==tissue)]<-tab$logFC
  }
}

if (length(remove)>0)
  
  fc.exp.mat<-fc.exp.mat[,-remove]

fc.exp.mat<-fc.exp.mat[rowSums(fc.exp.mat!=0)>0,]

labels<-read.table('/Users/karleg/Dimorph/labels.tsv',sep='\t',header=T)
    
for (i in (1:ncol(fc.exp.mat)))
  
  if (colnames(fc.exp.mat)[i] %in% labels$tissue)
    
    colnames(fc.exp.mat)[i]<-as.character(labels$X)[which(labels$tissue==colnames(fc.exp.mat)[i])]

colnames(fc.exp.mat)<-gsub('Brain - ','',colnames(fc.exp.mat))

colnames(fc.exp.mat)<-gsub(' \\(basal ganglia\\)','',colnames(fc.exp.mat))

colnames(fc.exp.mat)<-gsub(' \\(BA9\\)','',colnames(fc.exp.mat))

colnames(fc.exp.mat)<-gsub(' \\(BA24\\)','',colnames(fc.exp.mat))

colnames(fc.exp.mat)<-gsub(' \\(Suprapubic\\)','',colnames(fc.exp.mat))

colnames(fc.exp.mat)<-gsub(' \\(Lower leg\\)','',colnames(fc.exp.mat))

colnames(fc.exp.mat)<-gsub('Esophagus - Gastroesophageal Junction','Esophagus (GEJ)',colnames(fc.exp.mat))

dist.mat<-as.matrix(cor(fc.exp.mat))

colnames(dist.mat)<-colnames(fc.exp.mat)

rownames(dist.mat)<-colnames(fc.exp.mat)

library(pheatmap)

pheatmap(as.matrix(dist.mat))

save.image('/Users/karleg/Dimorph/RDATA/figure1b.RData')

