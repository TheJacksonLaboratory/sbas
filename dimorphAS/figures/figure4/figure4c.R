setwd('/Users/karleg/Dimorph/')

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

tissue.set<-tissue.sets[[16]]

comp<-tissue.set

event.type<-'se'

de.tab<-read.table(paste('/Users/karleg/Dimorph/gene_expression/DE_result_',comp,'.txt',sep=''))
      
as.tab<-read.table(paste('/Users/karleg/Dimorph/other/',paste(tissue.set,collapse='.'),event.type,'.txt',sep=''))
      
event.ensembl.map<-read.table(paste('/Users/karleg/Dimorph/fromGTF.',toupper(event.type),'.txt',sep=''),header=T)
      
event.ensembl.map<-event.ensembl.map[event.ensembl.map$ID %in% rownames(as.tab),]
      
rownames(event.ensembl.map)<-event.ensembl.map$ID
      
as.tab<-merge(as.tab,event.ensembl.map,by='row.names')
      
as.tab<-as.tab[as.tab$geneSymbol %in% rownames(de.tab),]
      
sig.events<-read.table(paste('/Users/karleg/Dimorph/other/',paste(tissue.set,collapse = '.'),event.type,'.txt',sep=''))
 
sig.events<-sig.events[sig.events$adj.P.Val<=0.05 & abs(sig.events$logFC)>=log2(1.5),]
      
de.tab<-de.tab[rownames(de.tab) %in% as.tab$geneSymbol,]
      
as.tab<-cbind(as.tab,rep(0,nrow(as.tab)))
      
colnames(as.tab)[ncol(as.tab)]<-'DE.logFC'
      
for (i in (1:nrow(as.tab)))
        
      as.tab[i,'DE.logFC']<-de.tab$logFC[which(rownames(de.tab)==as.tab$geneSymbol[i])]
      
      
x=as.tab$DE.logFC[as.tab$Row.names %in% rownames(sig.events)]
      
y=as.tab$logFC[as.tab$Row.names %in% rownames(sig.events)]
      
rbp.regulated<-read.table(paste('/Users/karleg/Dimorph/rbp_target/',paste(tissue.set,collapse='.'),'-',event.type,'rbp_events.txt',sep=''))
   
plot(x,y,main=paste(comp,event.type,sep=' '),xlab='logFC DE',ylab='logFC AS',col=ifelse(as.tab$Row.names[as.tab$Row.names %in% rownames(sig.events)] %in% rbp.regulated$V1,'blue','black'),pch=20)
      
y.upper<-y[y>0]
      
x.upper<-x[y>0]
      
y.lower<-y[y<0]
      
x.lower<-x[y<0]
      
fit.upper<-lm(y.upper ~ x.upper + I(x.upper^2))
      
lines(sort(x.upper), fitted(fit.upper)[order(x.upper)], col='red', type='l',cex=2) 
      
fit.lower<-lm(y.lower ~ x.lower + I(x.lower^2))
      
lines(sort(x.lower), fitted(fit.lower)[order(x.lower)], col='red', type='l',cex=2) 
      
save.image('/Users/karleg/Dimorph/RDATA/figure4c.RData')



# Removal all variables from workspace
rm(list=ls())
#Before running the following, use the Session menu to set working directory to source file location
load('figure4c.RData')
library(latex2exp)


npgBlue<- rgb(60/256,84/256,136/256,1)
npgRed <- rgb(220/256,0,0,0.5)
npgGreen <- rgb(0,160/256,135/256,1)
npgBrown <- rgb(126/256,97/256,72/256,1)

par(mar=c(5,5,1,1))

plot(x,y,main='',
     xlab=TeX("more in males $\\leftarrow$ Differential expression (log FC) $\\rightarrow$ more in females"),
     ylab=TeX("more in males $\\leftarrow$ Exon skipping (log FC) $\\rightarrow$ more in females"),
     cex.lab=1.6,
     cex.axis=1.5,
     col=ifelse(as.tab$Row.names[as.tab$Row.names %in% rownames(sig.events)] %in% rbp.regulated$V1,npgBlue,npgRed),
     pch=20)

y.upper<-y[y>0]
x.upper<-x[y>0]
y.lower<-y[y<0]
x.lower<-x[y<0]

fit.upper<-lm(y.upper ~ x.upper + I(x.upper^2))

lines(sort(x.upper), fitted(fit.upper)[order(x.upper)], col=npgBrown, type='l',lwd=1.5) 

fit.lower<-lm(y.lower ~ x.lower + I(x.lower^2))

lines(sort(x.lower), fitted(fit.lower)[order(x.lower)], col=npgBrown, type='l',lwd=1.5) 

abline(v=0,col="black",lty=2,lwd=2)
abline(h=0,col="black",lty=2,lwd=2)

text(1,3,"Type I",cex=2.5)
text(1,-3,"Type II",cex=2.5)

