

library(ggplot2)

inc.iso.counts<-read.csv('/Users/karleg/rMATs-GTEx/rmats_final.se.jc.ijc.txt',header=TRUE)

skip.iso.counts<-read.csv('/Users/karleg/rMATs-GTEx/rmats_final.se.jc.sjc.txt',header=TRUE)



s.tab<-read.csv('/Users/karleg/Downloads/SplicingCalculatorData/hg38_ensembl_splicing_calculator_output.tsv',sep='\t',header = T)

s.tab<-s.tab[s.tab$DONOR!='NaN' & s.tab$ACCEPTOR!='NaN',]

s.tab<-cbind(s.tab,paste(paste0('chr',s.tab$X.CHR),s.tab$BEGIN,s.tab$END,s.tab$STRAND,sep='-'))

colnames(s.tab)[ncol(s.tab)]<-'EXON_ID'

setwd('/Users/karleg/Dimorph/')

meta.data<-read.table('2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

tissue.set<-"Heart - Left Ventricle"

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

as.tab<-as.tab[as.tab$Row.names %in% rownames(sig.events),]


s.tab<-s.tab[!duplicated(s.tab$EXON_ID),] #remove this to count one score per transcript instead of per exon

as.tab1<-as.tab[as.tab$logFC<0 & as.tab$DE.logFC>0,]

group.1<-which((s.tab$END %in% as.tab1$exonEnd) & (s.tab$BEGIN %in% as.tab1$exonStart_0base) & (paste0('chr',s.tab$X.CHR) %in% as.tab1$chr) & (s.tab$STRAND %in% as.tab1$strand))


as.tab2<-as.tab[as.tab$logFC>0 & as.tab$DE.logFC>0,]

group.2<-which((s.tab$END %in% as.tab2$exonEnd) & (s.tab$BEGIN %in% as.tab2$exonStart_0base) & (paste0('chr',s.tab$X.CHR) %in% as.tab2$chr) & (s.tab$STRAND %in% as.tab2$strand))


scores.group.1<-s.tab$DONOR[group.1]+s.tab$ACCEPTOR[group.1]

scores.group.2<-s.tab$DONOR[group.2]+s.tab$ACCEPTOR[group.2]

t.test(scores.group.1,scores.group.2)

df<-data.frame(score=c(scores.group.1,scores.group.2),type=c(rep(1,length(scores.group.1)),rep(2,length(scores.group.2))))

pn<-ggplot(df, aes(x=type, y=score,group=type,fill=type))  + geom_boxplot(notch=TRUE) + theme_minimal()  
pn<- pn + theme(axis.text = element_text(size=0, hjust=0.5),axis.ticks.x=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y = element_text(face="bold",  size=22),
                axis.text.y = element_text(size=22),
                legend.position="none")   +ylab('Score')
my.pval = paste0("p==", t.test(scores.group.1,scores.group.2)$p.value)
pn <- pn + annotate("text", x=1.4, y=5, label= as.character(my.pval),parse=TRUE, fontface="bold", size=8) 
pn

#############  LOAD DATA FROM THE ABOVE
# Removal all variables from workspace
rm(list=ls())
#Before running the following, use the Session menu to set working directory to source file location
load('daniels_scoresHeartLV.RData')
library(ggplot2)
library(ggsci)


npgBlue<- rgb(60/256,84/256,136/256,1)
npgRed <- rgb(220/256,0,0,0.5)
npgColors <- c(npgBlue,npgRed)


df$type <- as.factor(df$type)  # colors dont work with df$type being numbers


pn<-ggplot(df, aes(x=type, y=score,fill=type, group=type))  + geom_boxplot(notch=TRUE, outlier.shape=NA) + theme_minimal()  
pn<-pn+ scale_y_continuous(limits = quantile(df$score, c(0.1, 0.9))) + guides(fill=FALSE)
pn<- pn + theme(axis.text = element_text(size=0, hjust=0.5),axis.ticks.x=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y = element_text(face="bold",  size=32),
                axis.text.x = element_text(size=36),
                axis.text.y = element_text(size=36),
                legend.position="none")   +ylab('Donor + Acceptor Score')
#pn <- pn +  scale_fill_npg() 
my.pval = paste0("p==", t.test(scores.group.1,scores.group.2)$p.value)
# p val is 2.05e-05
my.pval = "p==2.05 %*% 10 ^ -5"  ##
pn <- pn + annotate("text", x=1.5, y=21, label= as.character(my.pval),parse=TRUE, fontface="bold", size=12) 
pn
 
## get numbers
mean(df$score[df$type=="Type I"])
mean(df$score[df$type=="Type II"])
sum(df$type=="Type I")
sum(df$type=="Type II")


e.lengths.group1<-as.tab1$exonEnd-as.tab1$exonStart_0base

e.lengths.group2<-as.tab2$exonEnd-as.tab2$exonStart_0base

mean(e.lengths.group1)

mean(e.lengths.group2)

wilcox.test(e.lengths.group1,e.lengths.group2,paired = F)


df<-data.frame(len=c(e.lengths.group1,e.lengths.group2),type=c(rep(1,length(e.lengths.group1)),rep(2,length(e.lengths.group2))))
df$type <-  ifelse (df$type==1, "Type II", "Type I") 
df$type <- as.factor(df$type) 

pn<-ggplot(df, aes(x=type, y=len,fill=type, group=type))  + geom_boxplot(notch=TRUE, outlier.shape=NA) + theme_minimal()  
pn<-pn+ scale_y_continuous(limits = quantile(df$len, c(0.1, 0.9))) + guides(fill=FALSE)
pn<- pn + theme(axis.text = element_text(size=0, hjust=0.5),axis.ticks.x=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y = element_text(face="bold",  size=32),
                axis.text.x = element_text(size=36),
                axis.text.y = element_text(size=36),
                legend.position="none")   +ylab('Exon Length')
#pn <- pn +  scale_fill_npg() 
my.pval = paste0("p==", wilcox.test(e.lengths.group1,e.lengths.group2,paired = F)$p.value)
# p val is 3.9462838487075e-12
my.pval = "p==???3.95 %*% 10 ^ -12"  ##
pn <- pn + annotate("text", x=1.5, y=200, label= as.character(my.pval),parse=TRUE, fontface="bold", size=12) 
pn



