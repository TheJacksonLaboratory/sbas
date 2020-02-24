
library(ggplot2)

tissue='Heart - Left Ventricle'

events.table<-read.table(paste('/Users/karleg/Dimorph/other/',tissue,'se.txt',sep=''))

annot.table<-read.table('/Users/karleg/Dimorph/fromGTF.SE.txt',header=T)

rownames(annot.table)<-annot.table$ID

merged.table<-merge(events.table,annot.table,by='row.names')

all.genes<-read.table('/Users/karleg/Dimorph/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct',sep='\t',header=T,skip=2,colClasses = c(rep("character", 2), rep("NULL", 11688)))

all.genes<-all.genes[!duplicated(all.genes$Description),]

#which(!merged.table$geneSymbol %in% all.genes$Description)

merged.table<-merged.table[order(merged.table$adj.P.Val),]

norm.counts<-read.table(paste('/Users/karleg/Dimorph/matrices/',tissue,'se_splice_levels.txt',sep=''))

norm.counts<-norm.counts[rownames(norm.counts) %in% merged.table$Row.names,]

norm.counts<-norm.counts[order(match(rownames(norm.counts),merged.table$Row.names)),]

norm.inc.counts<-norm.counts[,1:(ncol(norm.counts)/2)]

meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

meta.data$Sample_Name_s<-gsub('-','\\.',meta.data$Sample_Name_s)

meta.data<-meta.data[meta.data$Sample_Name_s %in% colnames(norm.inc.counts),]

meta.data<-meta.data[order(match(meta.data$Sample_Name_s,colnames(norm.inc.counts))),]

sex<-meta.data$sex_s

norm.skip.counts<-norm.counts[,(ncol(norm.counts)/2+1):ncol(norm.counts)]

norm.skip.counts<-norm.skip.counts[merged.table$geneSymbol %in% all.genes$Description,]

norm.inc.counts<-norm.inc.counts[merged.table$geneSymbol %in% all.genes$Description,]

merged.table<-merged.table[merged.table$geneSymbol %in% all.genes$Description,]

sig.rows<-which(merged.table$adj.P.Val<=0.05 & abs(merged.table$logFC)>=log2(1.5))

sd.order<-order(apply(norm.inc.counts,1,sd)+apply(norm.skip.counts,1,sd))

sd.order<-sd.order[!(sd.order %in% which(merged.table$adj.P.Val<=0.05))]

if (length(sig.rows)>100)
{
  select.events<-c(sd.order[1:100],sig.rows[1:100])
  
}else{
  
  select.events<-c(sd.order[1:length(sig.rows)],sig.rows)
  
}
norm.inc.counts<-2^(norm.inc.counts)-0.5

norm.inc.counts[norm.inc.counts<0]=0

norm.inc.counts<-norm.inc.counts[select.events,]

norm.skip.counts<-2^(norm.skip.counts)-0.5

norm.skip.counts[norm.skip.counts<0]=0

norm.skip.counts<-norm.skip.counts[select.events,]

merged.table<-merged.table[select.events,]

sig.event<-9

df<-data.frame(event=c(rep('INC-male',sum(sex=='male')),rep('INC-female',sum(sex=='female')),rep('SKIP-male',sum(sex=='male')),rep('SKIP-female',sum(sex=='female'))),
               count=log2(c(as.numeric(norm.inc.counts[which(merged.table$adj.P.Val<=0.05)[sig.event],]),
                          as.numeric(norm.skip.counts[which(merged.table$adj.P.Val<=0.05)[sig.event],]))+0.5))

df$event <- factor(df$event,
                       levels = c('INC-male','SKIP-male','INC-female','SKIP-female'),ordered = TRUE)

p <- ggplot(df, aes(x=event, y=count,color=event)) + geom_boxplot()+ggtitle(paste(merged.table$geneSymbol[merged.table$adj.P.Val<=0.05][sig.event],tissue))+xlab('')+ylab('log2(count+0.5)')
p

save.image('/Users/karleg/Dimorph/RDATA/figureS2c.RData')

######################################################################################
######################################################################################
######################################################################################
# Removal all variables from workspace
rm(list=ls())
#Before running the following, use the Session menu to set working directory to source file location
load('figureS2c.RData')
library(ggplot2)
library(ggsci)

mylabs <-c('INC-M','SKIP-M','INC-F','SKIP-F')

#+ scale_x_discrete(labels= mylabs)


p <- ggplot(df, aes(x=event, y=count,fill=event)) + geom_boxplot(notch=TRUE)+ theme_minimal() + scale_fill_npg() 
p <- p + theme(axis.text = element_text(size=24, hjust=0.5),
               axis.title.x=element_blank(),
               axis.title.y = element_text(face="bold",  size=22),
               axis.text.y = element_text(size=36),
               legend.position="none")   +ylab(expression(log[2]~(counts+0.5))) + scale_x_discrete(labels= mylabs)
#p <- p +  ggtitle(paste(merged.table$geneSymbol[merged.table$adj.P.Val<=0.05][sig.event],tissue))+xlab('')+ylab('log2(count+0.5)')

p


