event.type='se'

vals.m<-rep(0,length(tissue.sets))

vals.f<-rep(0,length(tissue.sets))

num.genes<-rep(0,length(tissue.sets))

de.genes.to.fem<-rep(0,length(tissue.sets))

de.genes.to.mal<-rep(0,length(tissue.sets))

for (tissue.set in lapply(tissue.sets,paste,collapse='.'))
  
{
  for (tissue in unlist(strsplit(tissue.set,split = '\\.')))
  {
    if (!file.exists(paste0('/Users/karleg/Dimorph/gene_expression/all_genes/DE_result_',tissue,'.txt')))
      
      next
    
    de.tab<-read.table(paste0('/Users/karleg/Dimorph/gene_expression/all_genes/DE_result_',tissue,'.txt'))
    
  # mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host="www.ensembl.org")
    
  #  gene.loc<-getBM(attributes = c("hgnc_symbol","chromosome_name"),
                   # filters    = "hgnc_symbol",
                 #   values     = gsub('\\..*','',rownames(de.tab)), 
                  #  mart       = mart)
    
   # y.genes<-gene.loc$hgnc_symbol[gene.loc$chromosome_name=='Y']
    
   # de.tab<-de.tab[!(rownames(de.tab) %in% y.genes),]
    
    de.genes.to.fem[which(unlist(lapply(tissue.sets,paste,collapse='.'))==tissue.set)]<-de.genes.to.fem[which(unlist(lapply(tissue.sets,paste,collapse='.'))==tissue.set)]+sum(de.tab$logFC>=log2(1.5) & de.tab$adj.P.Val<=0.05)
    
    names(de.genes.to.fem)[which(unlist(lapply(tissue.sets,paste,collapse='.')==tissue.set))]<-tissue.set
    
    de.genes.to.mal[which(unlist(lapply(tissue.sets,paste,collapse='.'))==tissue.set)]<-de.genes.to.mal[which(unlist(lapply(tissue.sets,paste,collapse='.'))==tissue.set)]+sum(de.tab$logFC<=(-log2(1.5)) & de.tab$adj.P.Val<=0.05)
    
    names(de.genes.to.mal)[which(unlist(lapply(tissue.sets,paste,collapse='.')==tissue.set))]<-tissue.set
    
  }
  
  if (!file.exists(paste('/Users/karleg/Dimorph/other/',tissue.set,event.type,'.txt',sep='')))
  {
    vals.m[which(lapply(tissue.sets,paste,collapse='.')==tissue.set)]<-NA
    
    names(vals.m)[which(lapply(tissue.sets,paste,collapse='.')==tissue.set)]<-tissue.set
    
    vals.f[which(lapply(tissue.sets,paste,collapse='.')==tissue.set)]<-NA
    
    names(vals.f)[which(lapply(tissue.sets,paste,collapse='.')==tissue.set)]<-tissue.set
    
    num.genes[which(lapply(tissue.sets,paste,collapse='.')==tissue.set)]<-NA
    
    names(num.genes)[which(lapply(tissue.sets,paste,collapse='.')==tissue.set)]<-tissue.set
    
    de.genes.to.fem[which(lapply(tissue.sets,paste,collapse='.')==tissue.set)]<-NA
    
    names(de.genes.to.fem)[which(lapply(tissue.sets,paste,collapse='.')==tissue.set)]<-tissue.set
    
    de.genes.to.mal[which(lapply(tissue.sets,paste,collapse='.')==tissue.set)]<-NA
    
    names(de.genes.to.mal)[which(lapply(tissue.sets,paste,collapse='.')==tissue.set)]<-tissue.set
    
    next
  } 
  
  events.tab<-read.table(paste('/Users/karleg/Dimorph/other/',tissue.set,event.type,'.txt',sep=''))
  
  anno.tab<-read.table(paste0('/Users/karleg/Dimorph/fromGTF.',toupper(event.type),'.txt'),header = T)
  
  anno.tab<-anno.tab[anno.tab$ID %in% rownames(events.tab[abs(events.tab$logFC)>=log2(1.5) & events.tab$adj.P.Val<=0.05,]),]
  
  vals.m[which(unlist(lapply(tissue.sets,paste,collapse='.'))==tissue.set)]<-nrow(events.tab[events.tab$logFC>=log2(1.5) & events.tab$adj.P.Val<=0.05,])
  
  names(vals.m)[which(unlist(lapply(tissue.sets,paste,collapse='.'))==tissue.set)]<-tissue.set
  
  vals.f[which(unlist(lapply(tissue.sets,paste,collapse='.'))==tissue.set)]<-nrow(events.tab[events.tab$logFC<=(-log2(1.5)) & events.tab$adj.P.Val<=0.05,])
  
  names(vals.f)[which(unlist(lapply(tissue.sets,paste,collapse='.'))==tissue.set)]<-tissue.set
  
  num.genes[which(unlist(lapply(tissue.sets,paste,collapse='.'))==tissue.set)]<-ifelse(nrow(anno.tab)>0,length(unique(anno.tab$geneSymbol)),0)
  
  names(num.genes)[which(unlist(lapply(tissue.sets,paste,collapse='.')==tissue.set))]<-tissue.set
  
}


labels<-read.csv('/Users/karleg/Dimorph/labels.tsv',sep='\t')

library(ggplot2)

df <- data.frame(Log2DEGenes=log2(as.numeric(c(de.genes.to.mal,de.genes.to.fem))),Gender=c(rep('Male',length(de.genes.to.mal)),rep('Female',length(de.genes.to.fem))),Tissue=rep(labels[match(labels[,1],names(de.genes.to.fem)),2],2))

df<-df[!is.na(df$Log2DEGenes),]

df<-df[df$Log2DEGenes>0,]

ggplot(df, aes(fill=Gender, y=Log2DEGenes, x=Tissue)) + 
  geom_bar(position="dodge", stat="identity") + theme(axis.text.x=element_text(angle = -90, hjust = 0))+ scale_fill_manual("legend", values = c("Male" = "Blue", "Female" = "Red"))+ylab(bquote('Number of DE genes ('* ~Log[2]*')'))+ theme(axis.text.x=element_text(vjust=0.2))



df <- data.frame(Log2DEGenes=log2(as.numeric(c(vals.m,vals.f))),Gender=c(rep('Male',length(vals.m)),rep('Female',length(vals.f))),Tissue=rep(labels[match(labels[,1],names(vals.f)),2],2))

df<-df[!is.na(df$Log2DEGenes),]

df<-df[df$Log2DEGenes>0,]

ggplot(df, aes(fill=Gender, y=Log2DEGenes, x=Tissue)) + 
  geom_bar(position="dodge", stat="identity") + theme(axis.text.x=element_text(angle = -90, hjust = 0))+ scale_fill_manual("legend", values = c("Male" = "Blue", "Female" = "Red"))+ylab(bquote('Number of ASE events ('* ~Log[2]*')'))+ theme(axis.text.x=element_text(vjust=0.2))

save.image('/Users/karleg/Dimorph/RDATA/figureS1.RData')

#Before running the following, use the Session menu to set working directory to source file location
load('figureS1.RData')
library(ggplot2)
library(ggsci)

p<- ggplot(df, aes(fill=Gender, y=Log2DEGenes, x=Tissue)) + geom_bar(position="dodge", stat="identity") +  theme_minimal() + scale_fill_npg() 
p<-p+theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.2),
           axis.title.x = element_blank(),
           text = element_text(size=20),
           axis.text = element_text(size=20, hjust=0.5),
           axis.title.y = element_text(size=20),
           legend.title = element_blank(),
           legend.position="top")+ 
  #scale_fill_manual("legend", values = c("Male", "Female"))+
  ylab(bquote('Number of differentially expressed genes ('*log[2]*')'))
p
