library(limma)

library(multtest)

library(Biobase)

library(yarn)

library(edgeR)

args <- commandArgs(trailingOnly = TRUE)

inc.iso.counts.mem<-read.csv(paste('rmats_final.',args[2],'.jc.ijc.txt',sep=''),header=TRUE)

skip.iso.counts.mem<-read.csv(paste('rmats_final.',args[2],'.jc.sjc.txt',sep=''),header=TRUE)


meta.data<-read.table('2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

tissue.sets<-list(c("Brain - Cerebellum","Brain - Cerebellar Hemisphere"),
                  c("Brain - Caudate (basal ganglia)","Brain - Nucleus accumbens (basal ganglia)","Brain - Putamen (basal ganglia)"),
                  c("Brain - Cortex","Brain - Frontal Cortex (BA9)","Brain - Anterior cingulate cortex (BA24)"),
                  c("Brain - Amygdala","Brain - Hippocampus"),
                  c("Esophagus - Gastroesophageal Junction","Esophagus - Muscularis"),
                  c("Skin - Sun Exposed (Lower leg)","Skin - Not Sun Exposed (Suprapubic)")
                  )
for (singleton.set in unique(as.character(meta.data$body_site_s[!meta.data$body_site_s %in% unlist(tissue.sets)])))
  
  tissue.sets[length(tissue.sets)+1]<-singleton.set

for (tissue.set in tissue.sets[as.integer(args[1])])
{
  
meta.data<-read.table('2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

inc.iso.counts<-inc.iso.counts.mem

meta.data$Sample_Name_s<-gsub('-','\\.',meta.data$Sample_Name_s)

meta.data<-meta.data[meta.data$Sample_Name_s %in% colnames(inc.iso.counts),]

rownames(inc.iso.counts)<-inc.iso.counts$ID

inc.iso.counts<-inc.iso.counts[,-1]

inc.iso.counts<-inc.iso.counts[,order(match(colnames(inc.iso.counts),meta.data$Sample_Name_s))]

inc.iso.counts<-inc.iso.counts[,!grepl('11ILO',meta.data$Sample_Name_s) & (meta.data$body_site_s %in% tissue.set)]

skip.iso.counts<-skip.iso.counts.mem

rownames(skip.iso.counts)<-skip.iso.counts$ID

skip.iso.counts<-skip.iso.counts[,-1]

skip.iso.counts<-skip.iso.counts[,order(match(colnames(skip.iso.counts),meta.data$Sample_Name_s))]

skip.iso.counts<-skip.iso.counts[,!grepl('11ILO',meta.data$Sample_Name_s) & (meta.data$body_site_s %in% tissue.set)]

meta.data<-meta.data[!grepl('11ILO',meta.data$Sample_Name_s),]

meta.data<-meta.data[meta.data$body_site_s %in% tissue.set,]

if (sum(colnames(skip.iso.counts)==meta.data$Sample_Name_s)!=nrow(meta.data))
  
  quit(1)

if (sum(colnames(inc.iso.counts)==meta.data$Sample_Name_s)!=nrow(meta.data))
  
  quit(1)

exprs.mat <- matrix(data=cbind(as.matrix(skip.iso.counts),as.matrix(inc.iso.counts)),ncol=2*ncol(skip.iso.counts),nrow=nrow(skip.iso.counts),dimnames=list(1:nrow(skip.iso.counts),c(colnames(skip.iso.counts),paste(colnames(skip.iso.counts),'2',sep=''))))

keep.events<-rep(T,nrow(exprs.mat))

groups=rep(meta.data$sex_s,2)

groups<-groups[colSums(exprs.mat)>0]

exprs.mat<-exprs.mat[,colSums(exprs.mat)>0]

for (group in c('male','female'))

  keep.events<-keep.events & (rowSums(cpm(exprs.mat[,groups %in% group]) > 1) >= 0.25*min(table(groups)))

rm(exprs.mat)

inc.counts.male<-inc.iso.counts[,meta.data$sex_s=='male']

inc.counts.female<-inc.iso.counts[,meta.data$sex_s=='female']

skip.counts.male<-skip.iso.counts[,meta.data$sex_s=='male']

skip.counts.female<-skip.iso.counts[,meta.data$sex_s=='female']

inc.counts.mat<-cbind(inc.counts.male,inc.counts.female)

rownames(inc.counts.mat)<-rownames(inc.iso.counts)

skip.counts.mat<-cbind(skip.counts.male,skip.counts.female)

rownames(skip.counts.mat)<-rownames(skip.iso.counts)

counts.mat<-cbind(inc.counts.male,inc.counts.female,skip.counts.male,skip.counts.female)

rownames(counts.mat)<-rownames(inc.iso.counts)

counts.mat<-counts.mat[keep.events,]

isoform<-c(rep(1,ncol(inc.counts.male)+ncol(inc.counts.female)),rep(0,ncol(skip.counts.male)+ncol(skip.counts.female)))

sex<-c(rep(1,ncol(inc.counts.male)),rep(0,ncol(inc.counts.female)),rep(1,ncol(skip.counts.male)),rep(0,ncol(skip.counts.female)))

block=rep(c(1:ncol(inc.counts.male),(1+ncol(inc.counts.male)):(ncol(inc.counts.female)+ncol(inc.counts.male))),2)

keep.cols<-colSums(counts.mat)>0

isoform<-isoform[keep.cols]

sex<-sex[keep.cols]

block<-block[keep.cols]

counts.mat<-counts.mat[,keep.cols]

counts.mat<-DGEList(counts=counts.mat)

counts.mat<-calcNormFactors(counts.mat)

if (length(unique(sex))==1)
  
  next

design <- model.matrix(~sex+isoform+sex*isoform)

correlation =  duplicateCorrelation(voom(counts.mat)$E,design,block=block)$consensus.correlation

v <- voom(counts.mat, design,block = block, correlation = correlation)

fit <- lmFit(v, design,block = block, correlation = correlation)

fit <- eBayes(fit, robust=TRUE)

res=topTable(fit, coef='sex:isoform',number=nrow(counts.mat))

#hist(res$P.Value,main = paste(tissue.set,collapse=','))

write.table(res[res$adj.P.Val<=0.05 & abs(res$logFC)>=log2(1.5),],paste(paste(tissue.set,collapse='.'),args[2],'.txt',sep=''),row.names = T,col.names = T,quote = F)

write.table(fit$coefficients,paste(paste(tissue.set,collapse='.'),args[2],'_coefs.txt',sep=''),quote = F,row.names = F)


}
