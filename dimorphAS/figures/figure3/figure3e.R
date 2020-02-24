library(coda)
source("/Users/karleg/Downloads/DBDA2Eprograms/DBDA2E-utilities.R")

tissue.list<-c('Heart - Left Ventricle',
               'Breast - Mammary Tissue',
               'Brain - Cortex.Brain - Frontal Cortex (BA9).Brain - Anterior cingulate cortex (BA24)',
               'Adrenal Gland'
               ,'Adipose - Subcutaneous','Muscle - Skeletal','Thyroid','Cells - Transformed fibroblasts','Artery - Aorta','Skin - Sun Exposed (Lower leg).Skin - Not Sun Exposed (Suprapubic)')


all.genes<-read.table('/Users/karleg/Dimorph/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct',sep='\t',header=T,skip=2,colClasses = c(rep("character", 2), rep("NULL", 11688)))

all.genes<-all.genes[!duplicated(all.genes$Description),]

rbp.names<-unique(gsub('_.*','',list.files('/Users/karleg/Dimorph/RBP_PSSMs/')))

rbp.names<-rbp.names[rbp.names %in% all.genes$Description]

summary.tab<-matrix(ncol=7,nrow=0)

colnames(summary.tab)<-c('Event','Gene', 'Sig. RBPs','Sig. Gene Expression','Sig. Sex','Tissue','Dimorphic')

top.rbps<-rbp.names

df<-data.frame(coef=NULL,rbp=NULL,tissue=NULL)

for (tissue in c('Heart - Left Ventricle',
                 'Breast - Mammary Tissue',
                 'Brain - Cortex.Brain - Frontal Cortex (BA9).Brain - Anterior cingulate cortex (BA24)','Thyroid','Artery - Aorta'))
{
  
  load(paste('/Users/karleg/Dimorph/McmcMostVaryingMoreSigs_',tissue,'.Rdata',sep=''))
  
  hdi<-HPDinterval(mcmcCoda)  
  
  #diagMCMC( mcmcCoda , parName=c("beta3[101]") )  
  
  s <- summary(mcmcCoda)
  
  m <- s$statistics[,"Mean"]
  
  #hist(m[grepl('beta3',names(m))])
  
  #names(m)[grepl('beta2',names(m))][abs(m[grepl('beta2',names(m))])>1]
  events.tab<-read.table(paste('/Users/karleg/Dimorph/EventsTable_',tissue,'.txt',sep=''),header = T)
  
  beta2.mat<-matrix(nrow=nrow(events.tab),ncol=length(rbp.names))
  
  for (rbp in (1:length(rbp.names)))
    
    for (event in (1:nrow(events.tab)))
    {
      
      var.name<-paste0('beta2[',event,',',rbp,']')
      
      low<-hdi[[1]][rownames(hdi[[1]])==var.name][1]
      
      high<-hdi[[1]][rownames(hdi[[1]])==var.name][2]
      
      beta2.mat[event,rbp]<-m[grepl(paste0('beta2\\[',event,',',rbp,'\\]'),names(m))]
      
      if (low<0 && high>0)
        
        beta2.mat[event,rbp]<-0
      
    }
  colnames(beta2.mat)=rbp.names
  
  
  for (rbp in top.rbps)
    
    df<-rbind(df,cbind(beta2.mat[events.tab$adj.P.Val<=0.05,colnames(beta2.mat)==rbp],rep(rbp,nrow(beta2.mat)),rep(tissue,nrow(beta2.mat)))  )
  
  
  beta3.vec=matrix(nrow=nrow(events.tab),ncol=1)
  
  for (event in (1:nrow(events.tab)))
  {
    var.name<-paste0('beta3[',event,']')
    
    low<-hdi[[1]][rownames(hdi[[1]])==var.name][1]
    
    high<-hdi[[1]][rownames(hdi[[1]])==var.name][2]
    
    beta3.vec[event]<-m[grepl(paste0('beta3\\[',event,'\\]'),names(m))]
    
    if (low<0 && high>0)
      
      beta3.vec[event]<-0
  }
  
  #par(mfrow=c(1,2))
  
  #hist(beta3.vec[events.tab$adj.P.Val<=0.05],main=paste('Expression in Dimorphic',tissue),xlab='Posterior Mean')
  
  #hist(beta3.vec[events.tab$adj.P.Val>0.05],main=paste('Expression in Other',tissue),xlab='Posterior Mean')
  
  write.table(beta3.vec[events.tab$adj.P.Val<=0.05],paste0('/Users/karleg/Dimorph/',tissue,'_Expression_comparison_sig.txt'),quote = F)
  
  write.table(beta3.vec[events.tab$adj.P.Val>0.05],paste0('/Users/karleg/Dimorph/',tissue,'_Expression_comparison_other.txt'),quote = F)
  
  
  beta1.vec=matrix(nrow=nrow(events.tab),ncol=1)
  
  for (event in (1:nrow(events.tab)))
  {
    
    var.name<-paste0('beta1[',event,']')
    
    low<-hdi[[1]][rownames(hdi[[1]])==var.name][1]
    
    high<-hdi[[1]][rownames(hdi[[1]])==var.name][2]
    
    beta1.vec[event]<-m[grepl(paste0('beta1\\[',event,'\\]'),names(m))]
    
    if (low<0 && high>0)
      
      beta1.vec[event]<-0
  }
  
  #par(mfrow=c(1,2))
  
  #hist(beta1.vec[events.tab$adj.P.Val<=0.05],main=paste('Sex in Dimorphic',tissue),xlab='Posterior Mean')
  
  #hist(beta1.vec[events.tab$adj.P.Val>0.05],main=paste('Sex in Other',tissue),xlab='Posterior Mean')
  
  write.table(beta1.vec[events.tab$adj.P.Val<=0.05],paste0('/Users/karleg/Dimorph/',tissue,'_Sex_comparison_sig.txt'),quote = F)
  
  write.table(beta1.vec[events.tab$adj.P.Val>0.05],paste0('/Users/karleg/Dimorph/',tissue,'_Sex_comparison_other.txt'),quote = F)
  
  
  #read events table and compare dimorphic and non-dimorphic for each RBP
  
  if (sum(events.tab$adj.P.Val<=0.05)>1)
  {
    out.tab<-rbind(colMeans(beta2.mat[events.tab$adj.P.Val<=0.05,]),colMeans(beta2.mat[events.tab$adj.P.Val>0.05,]))
  }else{
    out.tab<-rbind(beta2.mat[events.tab$adj.P.Val<=0.05,],colMeans(beta2.mat[events.tab$adj.P.Val>0.05,]))
  }
  rownames(out.tab)<-c('Dimorphic','Other')
  
  write.table(out.tab,paste0('/Users/karleg/Dimorph/',tissue,'_RBP_comparison.txt'),quote = F)
  
  for (i in (1:nrow(events.tab)))
  {
    
    next.row<-rep('',6)
    #Event, Gene, Sig. RBBs,Gene Expression,Sex,Tissue,dimorphic
    
    next.row[1]<-as.character(events.tab$geneSymbol)[i]
    
    sig.rbps=''
    
    for (j in (1:length(rbp.names)))
    {
      if (beta2.mat[i,j]!=0)
      {
        
        if (next.row[2]!='')
          
          next.row[2]<-paste0(next.row[2],',')
        
        next.row[2]<-paste0(next.row[2],rbp.names[j],'(',round(beta2.mat[i,j],2),')')
      }
    }
    
    #if (beta3.vec[i]!=0)
    
    next.row[3]<-round(beta3.vec[i],2)
    
    # if (beta1.vec[i]!=0)
    
    next.row[4]<-round(beta1.vec[i],2)
    
    next.row[5]<-tissue
    
    next.row[6]<-ifelse(events.tab$adj.P.Val[i]<=0.05,'Yes','No')
    
    next.row=c(events.tab$Row.names[i],next.row)
    
    summary.tab<-rbind(summary.tab,next.row)
    
  }
  
  
}
#diagMCMC( mcmcCoda , parName=c("beta4[12,2]") ) 
write.table(summary.tab,'/Users/karleg/Dimorph/summary_hbm.txt',sep='\t',quote = F,row.names = F,col.names = T)

summary.tab<-summary.tab[summary.tab[,'Dimorphic']=='Yes',]

colnames(df)<-c('Coef','RBP','Tissue')

df$Coef<-as.numeric(as.character(df$Coef))

labels<-read.table('/Users/karleg/Dimorph/labels.tsv',sep='\t',header=T)

df$Tissue<-as.character(df$Tissue)

for (i in 1:nrow(df))
  
  if (df$Tissue[i] %in% labels$tissue)
    
    df$Tissue[i]<-as.character(labels$X[which(as.character(labels$tissue)==as.character(df$Tissue[i]))])

#df$Tissue[which(df$Tissue %in% labels$tissue)]<-labels$X[]
#df.plot<-df.plot[df.plot$Tissue %in% c('Heart (LV)','Cortex','Breast','Thyroid','Aorta'),]
library("ggsci")
library("ggplot2")
library("gridExtra")
library(grid)


df<-df[df$Coef!=0,]

#top.rbps<-names(sort(unlist(lapply(lapply(split(df$Coef[df$Coef!=0],df$RBP[df$Coef!=0]),abs),mean)),decreasing = T))[1:10]

sum.pos<-sort(unlist(lapply(lapply(split(df$Coef,df$RBP),'>',0),sum)),decreasing = T)

sum.neg<-sort(unlist(lapply(lapply(split(df$Coef,df$RBP),'<',0),sum)),decreasing = T)

hnrnp<-c("HNRNPA1",   "HNRNPA1L2", "HNRNPA2B1" ,"HNRNPC","HNRNPCL1" ,"HNRNPF","HNRNPH1","HNRNPH2","HNRNPK", 
         "HNRNPL","HNRNPM", "HNRNPU")

srsf<-c("SRSF1","SRSF10","SRSF2","SRSF7","SRSF9")

sum.pos<-sum.pos[order(names(sum.pos))]

sum.neg<-sum.neg[order(names(sum.neg))]

pos.rbps<-names(which(sum.pos/(sum.pos+sum.neg)>=0.75 & (sum.pos+sum.neg>quantile(sum.pos+sum.neg,0.2))))

neg.rbps<-names(which(sum.pos/(sum.pos+sum.neg)<=0.25 & (sum.pos+sum.neg>quantile(sum.pos+sum.neg,0.2))))

cs.rbps<-names(which(sum.pos/(sum.pos+sum.neg)>0.25 & sum.pos/(sum.pos+sum.neg)<0.75 & (sum.pos+sum.neg>quantile(sum.pos+sum.neg,0.2))))

df.counts<-data.frame(type=c('Skip','Inc','CS'),counts=c(length(pos.rbps),length(neg.rbps),length(cs.rbps)))

df$RBP<-as.character(df$RBP)

sort(sum.pos,decreasing = T)

pos.rbps

df1<-df[df$RBP %in% c("HNRNPH2", "HNRNPU","U2AF2"),]

pn1<-ggplot(df1,aes(factor(RBP),Coef)) + geom_violin(aes(fill='red')) + scale_fill_manual(values = '#4DBBD5FF') 
pn1 <- pn1 + theme_minimal() +  theme(text = element_text(size=20),
                                      axis.text = element_text(size=20, hjust=0.5),
                                      axis.title.x=element_blank(),
                                      axis.title.y = element_text(size=24),
                                      plot.title = element_text(hjust = 0.5),
                                      legend.position = "none") + ylab("") + labs(title="")+ylim(-2,2)+ geom_hline(yintercept=0)
#pn <- pn +  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)


sort(sum.neg,decreasing = T)

neg.rbps

sort(unlist(lapply(lapply(split(df$Coef,df$RBP),abs),sum)),decreasing = T)

df2<-df[df$RBP %in% c("YBX1", "SRSF9","MATR3"),]

pn2<-ggplot(df2,aes(factor(RBP),Coef)) +geom_violin(aes(fill='blue'))+ scale_fill_manual(values = '#00A087FF')
pn2 <- pn2 + theme_minimal() +  theme(text = element_text(size=20),
                                      axis.text = element_text(size=20, hjust=0.5),
                                      axis.title.x=element_blank(),
                                      axis.title.y = element_text(size=24),
                                      plot.title = element_text(hjust = 0.5),
                                      legend.position = "none") + ylab("") + labs(title="")+ylim(-1,1)+ geom_hline(yintercept=0)
#pn <- pn +  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)


df3<-df[df$RBP %in% c("HNRNPK", "HNRNPA1L2","SRSF7"),]

pn3<-ggplot(df3,aes(factor(RBP),Coef)) +geom_violin(aes(fill='green')) + scale_fill_manual(values = '#E64B35FF') 
pn3 <- pn3 + theme_minimal() +  theme(text = element_text(size=20),
                                      axis.text = element_text(size=20, hjust=0.5),
                                      axis.title.x=element_blank(),
                                      axis.title.y = element_text(size=24),
                                      plot.title = element_text(hjust = 0.5),
                                      legend.position = "none") + ylab("") + labs(title="")+ylim(-2,2)+ geom_hline(yintercept=0)
#pn <- pn +  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)


pn4<-ggplot(df.counts, aes(type, counts)) +  geom_bar(aes(fill = type), position = "dodge", stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ scale_fill_npg() + guides(fill=FALSE)+theme(axis.title.x=element_blank(),axis.title=element_text(size=18))



grid.arrange(pn1,pn2,pn3,pn4, nrow = 4,
             left = textGrob('Coefficient',gp = gpar(fontsize = 20), rot = 90,vjust=1))


rbp.tab<-matrix(ncol=3,nrow=0)

colnames(rbp.tab)<-c('RBP','Total','Mean')

rbp.coefs<-split(df$Coef,df$RBP)

for (rbp in names(sum.pos)[which(sum.pos+sum.neg>quantile(sum.pos+sum.neg,0.2))])
{
  
  dist.mean<-mean(unlist(rbp.coefs[rbp]))
  
  tot.int<-length(unlist(rbp.coefs[rbp]))
  
  rbp.tab<-rbind(rbp.tab,c(rbp,tot.int,dist.mean))
  
}

write.table(rbp.tab,'RBP_summary.txt',sep='\t',quote = F,row.names = F)

save.image('/Users/karleg/Dimorph/RDATA/figure3e.RData')

#############################################################################
#############################################################################
# Removal all variables from workspace
rm(list=ls())


load("figure3e.Rdata")
library(grid)
library(gridExtra)
library(ggplot2)
library(ggsci)



pn1<-ggplot(df1,aes(factor(RBP),Coef)) + geom_violin(aes(fill='red')) + scale_fill_manual(values = '#4DBBD5FF') 
pn1 <- pn1 + theme_minimal() +  theme(text = element_text(size=20),
                                      axis.text = element_text(size=20, hjust=0.5),
                                      axis.title.x=element_blank(),
                                      axis.title.y = element_text(size=24),
                                      plot.title = element_text(hjust = 0.5),
                                      legend.position = "none") + ylab("") + labs(title="")+ylim(-2,2)+ geom_hline(yintercept=0)
#pn <- pn +  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)


sort(sum.neg,decreasing = T)

neg.rbps

sort(unlist(lapply(lapply(split(df$Coef,df$RBP),abs),sum)),decreasing = T)

df2<-df[df$RBP %in% c("YBX1", "SRSF9","MATR3"),]

pn2<-ggplot(df2,aes(factor(RBP),Coef)) +geom_violin(aes(fill='blue'))+ scale_fill_manual(values = '#00A087FF')
pn2 <- pn2 + theme_minimal() +  theme(text = element_text(size=20),
                                      axis.text = element_text(size=20, hjust=0.5),
                                      axis.title.x=element_blank(),
                                      axis.title.y = element_text(size=24),
                                      plot.title = element_text(hjust = 0.5),
                                      legend.position = "none") + ylab("") + labs(title="")+ylim(-1,1)+ geom_hline(yintercept=0)
#pn <- pn +  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)


df3<-df[df$RBP %in% c("HNRNPK", "HNRNPA1L2","SRSF7"),]

pn3<-ggplot(df3,aes(factor(RBP),Coef)) +geom_violin(aes(fill='green')) + scale_fill_manual(values = '#E64B35FF') 
pn3 <- pn3 + theme_minimal() +  theme(text = element_text(size=20),
                                      axis.text = element_text(size=20, hjust=0.5),
                                      axis.title.x=element_blank(),
                                      axis.title.y = element_text(size=24),
                                      plot.title = element_text(hjust = 0.5),
                                      legend.position = "none") + ylab("") + labs(title="")+ylim(-2,2)+ geom_hline(yintercept=0)
#pn <- pn +  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)


pn4<-ggplot(df.counts, aes(type, counts)) +  geom_bar(aes(fill = type), position = "dodge", stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ scale_fill_npg() + guides(fill=FALSE)+theme(axis.title.x=element_blank(),axis.title=element_text(size=18))



grid.arrange(pn1,pn2,pn3,pn4, nrow = 4,
             left = textGrob('Coefficient',gp = gpar(fontsize = 20), rot = 90,vjust=1))

df1$facet = 1
df2$facet = 2
df3$facet = 3
df_c <- rbind(df1, df2, df3)
p <- ggplot(df_c,aes(x = factor(RBP), y = Coef, fill = factor(RBP))) + 
  geom_violin() + facet_wrap(~facet, nrow = 3, scales = "free") + 
  xlab("") + ylab("") + ylim(-2,2)+ geom_hline(yintercept=0) + scale_fill_npg() +
  theme_minimal() +  theme(text = element_text(size=20),
                                      axis.text = element_text(size=8, hjust=0.5),
                                      axis.title.x=element_blank(),
                                      axis.title.y = element_text(size=8),
                                      plot.title = element_text(hjust = 0.5),
                                      legend.position = "none", 
                           strip.text = element_blank())
  
pn4_new <- ggplot(df.counts, aes(type, counts)) +  
  geom_bar(fill = "#00008B",color="black", position = "dodge", stat="identity") + 
  geom_text(aes(x = type, y = counts + 10, label = paste(100 * round(counts/sum(counts), 3), "%", sep = "")), size = 3) +
  guides(fill=FALSE) +
  xlab("") + scale_y_continuous(breaks = c(0, 20, 40), limits = c(0, 60))+
  theme_minimal() +
  theme(
        axis.text = element_text(size = 8), 
        axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title = element_text(size = 10),
        axis.title.y = element_text(vjust = 5)
        )
pn4_new  
p_grid <- arrangeGrob(p,pn4_new, nrow = 2, heights = c(3, 1), 
             left = textGrob('Coefficient',gp = gpar(fontsize = 10), rot = 90,vjust=1, hjust = -0.2))

ggsave(file = "fig3e_test.pdf", plot = p_grid, height = 4, width = 3)
ggsave(file = "fig3e_test.png", plot = p_grid, height = 4, width = 3)

