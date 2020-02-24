
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

all.as.genes<-c()

all.de.genes<-c()

event.counts<-matrix(nrow=length(tissue.sets),ncol=5)

as.gene.counts<-matrix(nrow=length(tissue.sets),ncol=5)

colnames(event.counts)<-c('se','a3ss','a5ss','mxe','ri')

rownames(event.counts)<-unlist(lapply(tissue.sets,paste,collapse='.'))

colnames(as.gene.counts)<-c('se','a3ss','a5ss','mxe','ri')

rownames(as.gene.counts)<-unlist(lapply(tissue.sets,paste,collapse='.'))

for (tissue.set in tissue.sets)
{
  
  for (tissue in tissue.set)
  {
    
    if (!file.exists(paste0('/Users/karleg/Dimorph/gene_expression/all_genes/DE_result_',tissue,'.txt')))
      
      next
    
    de.tab<-read.table(paste0('/Users/karleg/Dimorph/gene_expression/all_genes/DE_result_',tissue,'.txt'))
    
    de.genes<-unique(rownames(de.tab)[abs(de.tab$logFC)>=log2(1.5) & de.tab$adj.P.Val<=0.05])
    
    all.de.genes<-unique(c(all.de.genes,de.genes))
  }
  
  for (event.type in colnames(event.counts))
  {
    
    if (!file.exists(paste0('/Users/karleg/Dimorph/other/',paste(tissue.set,collapse='.'),event.type,'.txt')))
      
      next
    
    as.tab<-read.table(paste0('/Users/karleg/Dimorph/other/',paste(tissue.set,collapse='.'),event.type,'.txt'))
    
    an.tab<-read.table(paste0('/Users/karleg/Dimorph/fromGTF.',toupper(event.type),'.txt'),header=T)
    
    event.counts[paste(tissue.set,collapse='.'),event.type]<-sum(as.tab$adj.P.Val<=0.05 & abs(as.tab$logFC)>=log2(1.5))
    
    merged.tab<-merge(as.tab,an.tab,by.x='row.names',by.y='ID')
    
    all.as.genes<-unique(c(all.as.genes,as.character(merged.tab$geneSymbol[merged.tab$adj.P.Val<=0.05 & abs(merged.tab$logFC)>=log2(1.5)])))
    
    as.gene.counts[paste(tissue.set,collapse='.'),event.type]<-length(unique(merged.tab$geneSymbol[merged.tab$adj.P.Val<=0.05 & abs(merged.tab$logFC)>=log2(1.5)]))
    
  }
  
}



all.genes<-c()

for (event.type in colnames(event.counts))
{
  an.tab<-read.table(paste0('/Users/karleg/Dimorph/fromGTF.',toupper(event.type),'.txt'),header=T)
  
  all.genes<-unique(c(all.genes,as.character(an.tab$geneSymbol)))
  
}

de.univ<-sum(all.de.genes %in% all.genes)

de.as<-length(intersect(all.as.genes,all.de.genes))

total.univ<-length(all.genes)

total.as<-length(all.as.genes)

exp(phyper(q=de.as, m=de.univ, n=total.univ-de.univ, k=total.as,log.p = T,lower.tail = F))

library(VennDiagram)

venn.plot <- draw.triple.venn(
  area1 = de.univ,
  area2 = total.as,
  area3 = total.univ,
  n12 = de.as,
  n23 = total.as,
  n13 = de.univ,
  n123 = de.as,
  category = c("DE", "sex-biased AS", "ALL AS"),
  fill = c("blue", "red", "green"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("blue", "red", "gold")
)
grid.draw(venn.plot)

save.image('/Users/karleg/Dimorph/RDATA/venn_diagram.RData')

###########################################################
##########################################################
load('venn_diagram.RData')
library(grid)
library(VennDiagram)


venn.plot <- draw.triple.venn(
  area1 = de.univ,
  area2 = total.as,
  area3 = total.univ,
  n12 = de.as,
  n23 = total.as,
  n13 = de.univ,
  n123 = de.as,
  category = c("DE", "sex-biased AS", "ALL AS"),
  fill = c("blue", "red", "green"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("blue", "red", "gold")
)
grid.draw(venn.plot)
