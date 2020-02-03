


event.type='se'

meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

tissue.sets<-list(c("Brain - Cerebellum","Brain - Cerebellar Hemisphere"),
                  c("Brain - Caudate (basal ganglia)","Brain - Nucleus accumbens (basal ganglia)","Brain - Putamen (basal ganglia)"),
                  c("Brain - Cortex","Brain - Frontal Cortex (BA9)","Brain - Anterior cingulate cortex (BA24)"),
                  c("Brain - Amygdala","Brain - Hippocampus"),
                  c("Esophagus - Gastroesophageal Junction","Esophagus - Muscularis"),
                  c("Skin - Sun Exposed (Lower leg)","Skin - Not Sun Exposed (Suprapubic)")
)
for (singleton.set in unique(as.character(meta.data$body_site_s[!meta.data$body_site_s %in% unlist(tissue.sets)])))
  
  tissue.sets[length(tissue.sets)+1]<-singleton.set


an.tab<-read.table(paste0('/Users/karleg/Dimorph/fromGTF.',toupper(event.type),'.txt'),header=T)

gene.names<-as.character(an.tab$geneSymbol)

fc.mat<-matrix(rep(0), nrow(an.tab)*length(tissue.sets),nrow = nrow(an.tab),ncol=length(tissue.sets))

colnames(fc.mat)<-lapply(tissue.sets,paste,collapse='.')

remove=c()

for (tissue.set in tissue.sets)
{
  
  if (!file.exists(paste0('/Users/karleg/Dimorph/other/',paste(tissue.set,collapse='.'),tolower(event.type),'.txt')))
  {
    remove<-c(remove,which(colnames(fc.mat)==paste(tissue.set,collapse='.')))
    
  }else{
    
    tab<-read.table(paste0('/Users/karleg/Dimorph/other/',paste(tissue.set,collapse='.'),tolower(event.type),'.txt'))
    
    fc.mat[as.integer(rownames(tab)),which(colnames(fc.mat)==paste(tissue.set,collapse='.'))]<-tab$logFC
  }
}

fc.mat<-fc.mat[,-remove]

gene.names<-gene.names[rowSums(fc.mat!=0)>0]

fc.mat<-fc.mat[rowSums(fc.mat!=0)>0,]

fc.mat<-as.data.frame(fc.mat)

fc.mat<-do.call(rbind,lapply(split(fc.mat,gene.names),colMeans))

fc.mat<-normalizeQuantiles(fc.mat)

dist.mat<-as.matrix(cor(fc.mat))

colnames(dist.mat)<-colnames(fc.mat)

rownames(dist.mat)<-colnames(fc.mat)

library(pheatmap)

pheatmap(as.matrix(dist.mat))



