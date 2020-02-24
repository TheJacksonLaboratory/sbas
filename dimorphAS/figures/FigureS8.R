library(multtest)

library(pheatmap)
#  number of binding sites in significant events for each rbp for each tissue
#  normalized mean expression values for each rbp over the tissues (normalize by dividing by the max expression level across all tissues)
#can be visualized using a stacked bar plot

setwd('/Users/karleg/Dimorph/')

meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

event.type='se'

tissue.sets<-list(c("Brain - Cerebellum","Brain - Cerebellar Hemisphere"),
                  c("Brain - Caudate (basal ganglia)","Brain - Nucleus accumbens (basal ganglia)","Brain - Putamen (basal ganglia)"),
                  c("Brain - Cortex","Brain - Frontal Cortex (BA9)","Brain - Anterior cingulate cortex (BA24)"),
                  c("Brain - Amygdala","Brain - Hippocampus"),
                  c("Esophagus - Gastroesophageal Junction","Esophagus - Muscularis"),
                  c("Skin - Sun Exposed (Lower leg)","Skin - Not Sun Exposed (Suprapubic)")
)
for (singleton.set in unique(as.character(meta.data$body_site_s[!meta.data$body_site_s %in% unlist(tissue.sets)])))
  
  tissue.sets[length(tissue.sets)+1]<-singleton.set


clean.sets<-list()

for (i in (1:length(tissue.sets)))
{
  if (!file.exists(paste('/Users/karleg/Dimorph/other/',paste(tissue.sets[[i]],collapse='.'),event.type,'.txt',sep='')))
  { 
    remove.tissue<-T
  }else{
    sig.events<-read.table(paste('/Users/karleg/Dimorph/other/',paste(tissue.sets[[i]],collapse = '.'),event.type,'.txt',sep=''),header=T)
    
    sig.events<-sig.events[sig.events$adj.P.Val<=0.05 & abs(sig.events$logFC)>=log2(1.5),]
    
    remove.tissue<-nrow(sig.events)<=5
  }
  
  if (!remove.tissue)
    
    clean.sets[[length(clean.sets)+1]]<-tissue.sets[[i]]
}

tissue.sets<-clean.sets

pwm.files<-list.files('RBP_PSSMs/')  #files with the PSSMs for the different RBP

factor.names<-unique(gsub('_.*','',pwm.files))

if (!file.exists('rbp_expression.txt'))
  
  system(paste('cat GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct | grep -e Name -e ',paste(paste('[[:space:]]',factor.names,'[[:space:]]',sep=''),sep='',collapse = ' -e '),' > rbp_expression.txt',sep=''))

rbp.exp<-read.table('rbp_expression.txt',header=TRUE)

n.rbps<-nrow(rbp.exp) #number of RBP in the expression data

rbp.exp<-t(rbp.exp)

colnames(rbp.exp)<-rbp.exp[2,]

rbp.exp<-rbp.exp[-c(1,2),]

expressed.rbps<-colnames(rbp.exp)

rownames(meta.data)<-gsub('-','\\.',meta.data$Sample_Name_s)

rbp.exp<-merge(x=rbp.exp,y=meta.data,by='row.names')

inc.iso.counts<-read.csv(paste('/Users/karleg/rMATs-GTEx/rmats_final.se.jc.ijc.txt',sep=''),header=TRUE,nrows = 1)

rbp.exp<-rbp.exp[rbp.exp$Row.names %in% colnames(inc.iso.counts),]

rbp.exp<-rbp.exp[!grepl('11ILO',rbp.exp$Row.names),]



skip.iso.counts.mem<-read.csv(paste('/Users/karleg/rMATs-GTEx/rmats_final.',event.type,'.jc.sjc.txt',sep=''),header=TRUE)


event.type='se'

tissue.set<-"Breast - Mammary Tissue"

sig.events<-read.table(paste('/Users/karleg/Dimorph/other/',paste(tissue.set,collapse='.'),event.type,'.txt',sep=''),header=T)

sig.events<-sig.events[sig.events$adj.P.Val<=0.05 & abs(sig.events$logFC)>=log2(1.5),]

skip.iso.counts<-skip.iso.counts.mem

rownames(skip.iso.counts)<-skip.iso.counts$ID

skip.iso.counts<-skip.iso.counts[,-1]

colnames(skip.iso.counts)<-gsub('\\.','-',colnames(skip.iso.counts))

skip.iso.counts<-skip.iso.counts[,colnames(skip.iso.counts) %in% rbp.exp$Sample_Name_s]

skip.iso.counts<-skip.iso.counts[,order(match(colnames(skip.iso.counts),rbp.exp$Sample_Name_s))]

skip.iso.counts<-skip.iso.counts[,(rbp.exp$body_site_s %in% tissue.set)]

snrnp70.mat<-skip.iso.counts[rownames(skip.iso.counts) %in% rownames(sig.events),]

factor.n<-'FUS'

design<-as.matrix(as.numeric(rbp.exp[,colnames(rbp.exp)==factor.n]),ncol=1)

design<-as.matrix(design[rbp.exp$body_site_s %in% tissue.set,],ncol=1)

res=c()

if (nrow(snrnp70.mat)==0)
  
  next

direction=c()

for (i in (1:nrow(snrnp70.mat)))
{
  
  fit<-lm(as.numeric(snrnp70.mat[i,])~design)
  
  res<-c(res,summary(fit)$coefficients[2,4])
  
  direction<-c(direction,sign(summary(fit)$coefficients[2,1]))
  
}

c.res<-sum(res<=0.05)

if (length(res)>1)
  
  c.res<-mt.rawp2adjp(res,proc='BH')[[1]][,2][order(mt.rawp2adjp(res,proc='BH')$index)]

which.dir=-1

fit <- lm(colMeans(snrnp70.mat[!is.na(c.res) & c.res<=0.05 & direction==which.dir,])~design)

plot(design,colMeans(snrnp70.mat[!is.na(c.res) & c.res<=0.05 & direction==which.dir,]),ylab='skip event mean counts',xlab=paste(factor.n,'expression'),main=paste(factor.n,ifelse(which.dir==1,'Positive','Negative'),paste(tissue.set,collapse='.')))

rmse <- round(sqrt(mean(resid(fit)^2)), 2)

coefs <- coef(fit)

b0 <- round(coefs[1], 4)

b1 <- round(coefs[2],4)

r2 <- round(summary(fit)$r.squared, 2)

eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ r^2 == .(r2) * "," ~~ RMSE == .(rmse) ~~ Fvalp==.(1-pf(summary(fit)$fstatistic[1],df1=summary(fit)$fstatistic[2],df2=summary(fit)$fstatistic[3])))

abline(fit,col='red')

mtext(side=1,eqn) 



target.event<-10565

fit <- lm(colMeans(snrnp70.mat[rownames(snrnp70.mat)==target.event,])~design)

plot(design,colMeans(snrnp70.mat[rownames(snrnp70.mat)==target.event,]),ylab='skip event mean counts',xlab=paste(factor.n,'expression'),main=paste(factor.n,target.event,paste(tissue.set,collapse='.')))

rmse <- round(sqrt(mean(resid(fit)^2)), 2)

coefs <- coef(fit)

b0 <- round(coefs[1], 4)

b1 <- round(coefs[2],4)

r2 <- round(summary(fit)$r.squared, 2)

eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ r^2 == .(r2) * "," ~~ RMSE == .(rmse) ~~ Fvalp==.(1-pf(summary(fit)$fstatistic[1],df1=summary(fit)$fstatistic[2],df2=summary(fit)$fstatistic[3])))

abline(fit,col='red')

mtext(side=1,eqn) 


for (factor.n in c('HNRNPA1','FUS','SNRNP70','YBX1'))
  
  for (target.event in c('8160','10565','21944','22731'))
  {
    
    design<-as.matrix(as.numeric(rbp.exp[,colnames(rbp.exp)==factor.n]),ncol=1)
    
    design<-as.matrix(design[rbp.exp$body_site_s %in% tissue.set,],ncol=1)
    
    fit <- lm(colMeans(snrnp70.mat[rownames(snrnp70.mat)==target.event,])~design)
    
    plot(design,colMeans(snrnp70.mat[rownames(snrnp70.mat)==target.event,]),ylab='skip event mean counts',xlab=paste(factor.n,'expression'),main=paste(factor.n,target.event,paste(tissue.set,collapse='.')))
    
    rmse <- round(sqrt(mean(resid(fit)^2)), 2)
    
    coefs <- coef(fit)
    
    b0 <- round(coefs[1], 4)
    
    b1 <- round(coefs[2],4)
    
    r2 <- round(summary(fit)$r.squared, 2)
    
    eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ r^2 == .(r2) * "," ~~ RMSE == .(rmse) ~~ Fvalp==.(1-pf(summary(fit)$fstatistic[1],df1=summary(fit)$fstatistic[2],df2=summary(fit)$fstatistic[3])))
    
    abline(fit,col='red')
    
    mtext(side=1,eqn) 
    
  }
      
   

for (which.dir in c(1,-1))
  for (factor.n in c('HNRNPA1','FUS','SNRNP70','YBX1'))
  {
    
    design<-as.matrix(as.numeric(rbp.exp[,colnames(rbp.exp)==factor.n]),ncol=1)
    
    design<-as.matrix(design[rbp.exp$body_site_s %in% tissue.set,],ncol=1)
    
    res=c()
    
    direction=c()
    
    for (i in (1:nrow(snrnp70.mat)))
    {
      
      fit<-lm(as.numeric(snrnp70.mat[i,])~design)
      
      res<-c(res,summary(fit)$coefficients[2,4])
      
      direction<-c(direction,sign(summary(fit)$coefficients[2,1]))
      
    }
    
    c.res<-sum(res<=0.05)
    
    if (length(res)>1)
      
      c.res<-mt.rawp2adjp(res,proc='BH')[[1]][,2][order(mt.rawp2adjp(res,proc='BH')$index)]
    
    fit <- lm(colMeans(snrnp70.mat[!is.na(c.res) & c.res<=0.05 & direction==which.dir,])~design)
    
    plot(design,colMeans(snrnp70.mat[!is.na(c.res) & c.res<=0.05 & direction==which.dir,]),ylab='skip event mean counts',xlab=paste(factor.n,'expression'),main=paste(factor.n,ifelse(which.dir==1,'Positive','Negative'),paste(tissue.set,collapse='.')))
    
    rmse <- round(sqrt(mean(resid(fit)^2)), 2)
    
    coefs <- coef(fit)
    
    b0 <- round(coefs[1], 4)
    
    b1 <- round(coefs[2],4)
    
    r2 <- round(summary(fit)$r.squared, 2)
    
    eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ r^2 == .(r2) * "," ~~ RMSE == .(rmse) ~~ Fvalp==.(1-pf(summary(fit)$fstatistic[1],df1=summary(fit)$fstatistic[2],df2=summary(fit)$fstatistic[3])))
    
    abline(fit,col='red')
    
    mtext(side=1,eqn) 
    
    
  
  }
  
  
save.image('/Users/karleg/Dimorph/RDATA/figureS8.RData')  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  


