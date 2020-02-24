library(Biobase)

library(yarn)

#Merge if for all event types the diatance is >=0.95
meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

screen.list<-list(grepl('Skin',meta.data$body_site_s) | grepl('fibro',meta.data$body_site_s),grepl('Colon',meta.data$body_site_s),grepl('Artery',meta.data$body_site_s),grepl('Esophagus',meta.data$body_site_s),grepl('Adipose',meta.data$body_site_s))

for (screen in screen.list)
{
tissue.sets<-as.character(unique(meta.data$body_site_s[screen]))  #Colon,Skin and fibroblasts,Artery, Esophagus,Brain

all.events.mat<-matrix(nrow=length(tissue.sets),ncol=length(tissue.sets),data=rep(1,length(tissue.sets)^2))

rownames(all.events.mat)<-1:nrow(all.events.mat)

colnames(all.events.mat)<-1:ncol(all.events.mat)


for (event.type in c('se','mxe','ri','a5ss','a3ss'))
{
inc.iso.counts<-read.csv(paste('/Users/karleg/rMATs-GTEx/rmats_final.',event.type,'.jc.ijc.txt',sep=''),header=TRUE)

rownames(inc.iso.counts)<-inc.iso.counts$ID

inc.iso.counts<-inc.iso.counts[,-1]

meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

meta.data$Sample_Name_s<-gsub('-','\\.',meta.data$Sample_Name_s)

meta.data<-meta.data[meta.data$Sample_Name_s %in% colnames(inc.iso.counts),]

skip.iso.counts<-read.csv(paste('/Users/karleg/rMATs-GTEx/rmats_final.',event.type,'.jc.sjc.txt',sep=''),header=TRUE)

rownames(skip.iso.counts)<-skip.iso.counts$ID

skip.iso.counts<-skip.iso.counts[,-1]

diff.mat<-matrix(nrow=length(tissue.sets),ncol=length(tissue.sets))

rownames(diff.mat)=1:nrow(diff.mat)

colnames(diff.mat)=1:ncol(diff.mat)

for (i in (1:ncol(combn(length(tissue.sets),2))))
{
  
tissue.set<-as.character(tissue.sets[combn(length(tissue.sets),2)[,i]] )
  
rownames(diff.mat)[combn(length(tissue.sets),2)[1,i]]=tissue.set[1]

rownames(diff.mat)[combn(length(tissue.sets),2)[2,i]]=tissue.set[2]

colnames(diff.mat)[combn(length(tissue.sets),2)[1,i]]=tissue.set[1]

colnames(diff.mat)[combn(length(tissue.sets),2)[2,i]]=tissue.set[2]
  
  
meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)
  
meta.data$Sample_Name_s<-gsub('-','\\.',meta.data$Sample_Name_s)

meta.data<-meta.data[meta.data$Sample_Name_s %in% colnames(inc.iso.counts),]

inc.iso.counts<-inc.iso.counts[order(match(colnames(inc.iso.counts),meta.data$Sample_Name_s))]

skip.iso.counts<-skip.iso.counts[order(match(colnames(skip.iso.counts),meta.data$Sample_Name_s))]

if (sum(colnames(inc.iso.counts)==meta.data$Sample_Name_s)!=nrow(meta.data) || sum(colnames(skip.iso.counts)==meta.data$Sample_Name_s)!=nrow(meta.data))
  
  quit(1)



tissue.skip<-skip.iso.counts[,meta.data$body_site_s %in% tissue.set]

tissue.inc<-inc.iso.counts[,meta.data$body_site_s %in% tissue.set]

exprs.mat <- rbind(as.matrix(tissue.skip),as.matrix(tissue.inc))

colnames(exprs.mat)<-colnames(tissue.skip)

rownames(exprs.mat)<-1:nrow(exprs.mat)

pData <- data.frame(row.names=colnames(tissue.skip),sex=meta.data$sex_s[meta.data$body_site_s %in% tissue.set])

phenoData <- AnnotatedDataFrame(data=pData)

obj <- ExpressionSet(exprs.mat, phenoData=phenoData)

res=checkTissuesToMerge(obj,as.character(meta.data$histological_type_s[meta.data$body_site_s %in% tissue.set]),as.factor(as.character(meta.data$body_site_s[meta.data$body_site_s %in% tissue.set])),legendPos='bottom')

dist.1<-mean(dist(res[[1]][meta.data$body_site_s[meta.data$body_site_s %in% tissue.set]==tissue.set[1],]))

dist.2<-mean(dist(res[[1]][meta.data$body_site_s[meta.data$body_site_s %in% tissue.set]==tissue.set[2],]))

dist.all<-mean(dist(res[[1]]))

diff.mat[combn(length(tissue.sets),2)[1,i],combn(length(tissue.sets),2)[2,i]]=(dist.1+dist.2)/(2*dist.all)

diff.mat[combn(length(tissue.sets),2)[2,i],combn(length(tissue.sets),2)[1,i]]=(dist.1+dist.2)/(2*dist.all)

res=checkTissuesToMerge(obj[,meta.data$body_site_s[meta.data$body_site_s %in% tissue.set] %in% tissue.set[1]],as.character(meta.data$histological_type_s[meta.data$body_site_s %in% tissue.set[1]]),as.factor(as.character(meta.data$body_site_s[meta.data$body_site_s %in% tissue.set[1]])),legendPos='bottom')

dist.1<-mean(dist(res[[1]]))

rand.res=c()

for (j in 1:1000)
{
rand.samp<-sample(c(TRUE,FALSE),nrow(res[[1]]),replace = T)

rdist.1<-dist(res[[1]][ rand.samp,])

rdist.2<-rdist.1<-dist(res[[1]][!rand.samp,])

rand.res=c(rand.res,(mean(rdist.1)+mean(rdist.2))/(2*dist.1))
}
#hist(rand.res)

}

for (i in (1:nrow(diff.mat)))

  diff.mat[i,i]=1


diff.mat[diff.mat>=0.95]<-1

diff.mat[diff.mat<0.95]<-0

all.events.mat<-all.events.mat*diff.mat

#library(pheatmap)

#pheatmap(diff.mat)
meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)


}

rownames(all.events.mat)<-rownames(diff.mat)

colnames(all.events.mat)<-colnames(diff.mat)

jpeg(file=paste('/Users/karleg/Dimorph/figures/',paste(tissue.set,collapse='.'),sep=''))

pheatmap(all.events.mat)

meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

dev.off()
}
#Sex annotation

meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

meta.data$Sample_Name_s<-gsub('-','\\.',meta.data$Sample_Name_s)

meta.data<-meta.data[meta.data$Sample_Name_s %in% colnames(inc.iso.counts),]

tissue.sets<-unique(meta.data$body_site_s)

for (tissue.set in tissue.sets)
{

meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)
  
meta.data$Sample_Name_s<-gsub('-','\\.',meta.data$Sample_Name_s)
  
meta.data<-meta.data[meta.data$Sample_Name_s %in% colnames(inc.iso.counts),]
  
tissue.skip<-skip.iso.counts[,(meta.data$body_site_s %in% tissue.set) & !grepl('11ILO',meta.data$Sample_Name_s)]

tissue.inc<-inc.iso.counts[,(meta.data$body_site_s %in% tissue.set) & !grepl('11ILO',meta.data$Sample_Name_s)]

meta.data<-meta.data[!grepl('11ILO',meta.data$Sample_Name_s),]

exprs.mat <- rbind(as.matrix(tissue.skip),as.matrix(tissue.inc))

colnames(exprs.mat)<-colnames(tissue.skip)

rownames(exprs.mat)<-1:nrow(exprs.mat)

pData <- data.frame(row.names=colnames(tissue.skip),sex=meta.data$sex_s[meta.data$body_site_s %in% tissue.set])

phenoData <- AnnotatedDataFrame(data=pData)

fData<-data.frame(chr=rep(read.table('/Users/karleg/Dimorph/fromGTF.SE.txt',header=T)$chr,2),row.names =rownames(exprs.mat) )

featureData<-AnnotatedDataFrame(data=fData)

obj <- ExpressionSet(exprs.mat, phenoData=phenoData,featureData = featureData)

res=checkMisAnnotation(obj,'sex',controlGenes = 'chrY',columnID = 'chr',main=tissue.set)

text(res[,1],res[,2],labels=lapply(strsplit(colnames(tissue.skip),'\\.'),'[',2))
}
