tab<-read.table('/Users/karleg/Downloads/nature24265-s3/Suppl.Table.1.txt',sep='\t',header=T)

table(tab$Combined.XCI.status[!is.na(tab$Gene.name)])

res.exp<-matrix(ncol=4,nrow=0)

colnames(res.exp)<-c('Tissue','Escape-Biased','Inactive-Biased','Variable-Biased')

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

for (tissue in unique(meta.data$body_site_s))
{
  
  if (!file.exists(paste0('/Users/karleg/Dimorph/gene_expression/all_genes/DE_result_',tissue,'.txt')))
    
    next
  
  de.tab<-read.table(paste0('/Users/karleg/Dimorph/gene_expression/all_genes/DE_result_',tissue,'.txt'))
  
  de.genes<-unique(rownames(de.tab)[abs(de.tab$logFC)>=log2(1.5) & de.tab$adj.P.Val<=0.05])
  
  escape.biased<-sum(unique(as.character(tab$Gene)[tab$Combined.XCI.status=='escape']) %in% unique(de.genes))
  
  inactive.biased<-sum(unique(as.character(tab$Gene)[tab$Combined.XCI.status=='inactive']) %in% unique(de.genes))
  
  variable.biased<-sum(unique(as.character(tab$Gene)[tab$Combined.XCI.status=='variable']) %in% unique(de.genes))
  
 # escape.biased<-escape.biased/table(tab$Combined.XCI.status[tab$Gene.name %in% rownames(de.tab)])['escape']
  
 # inactive.biased<-inactive.biased/table(tab$Combined.XCI.status[tab$Gene.name %in% rownames(de.tab)])['inactive']
  
 # variable.biased<-inactive.biased/table(tab$Combined.XCI.status[tab$Gene.name %in% rownames(de.tab)])['variable']
  
  res.exp<-rbind(res.exp,c(tissue,escape.biased,inactive.biased,variable.biased))
  
}

res.exp<-rbind(res.exp,c('Sum',colSums(apply(res.exp[,2:4],2,as.numeric))))

#write.table(res.exp,'/Users/karleg/Dimorph/XiDEgenes.txt',sep='\t',row.names = F,quote = F)

res.as<-matrix(ncol=4,nrow=0)

colnames(res.as)<-c('Tissue','Escape-Biased','Inactive-Biased','Variable-Biased')

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

event.type<-'se'

for (tissue.set in tissue.sets)
{
 
  if (!file.exists(paste0('/Users/karleg/Dimorph/other/',paste(tissue.set,collapse='.'),event.type,'.txt')))
    
    next
  
  as.tab<-read.table(paste0('/Users/karleg/Dimorph/other/',paste(tissue.set,collapse='.'),event.type,'.txt'))
  
  an.tab<-read.table(paste0('/Users/karleg/Dimorph/fromGTF.',toupper(event.type),'.txt'),header=T)
  
  merged.tab<-merge(as.tab,an.tab,by.x='row.names',by.y='ID')
  
  as.genes<-as.character(merged.tab$geneSymbol[abs(merged.tab$logFC)>=log2(1.5) & merged.tab$adj.P.Val<=0.05])
  
  escape.biased<-sum(unique(as.character(tab$Gene)[tab$Combined.XCI.status=='escape']) %in% unique(as.genes))
  
  inactive.biased<-sum(unique(as.character(tab$Gene)[tab$Combined.XCI.status=='inactive']) %in% unique(as.genes))
  
  variable.biased<-sum(unique(as.character(tab$Gene)[tab$Combined.XCI.status=='variable']) %in% unique(as.genes))

 # escape.biased<-escape.biased/table(tab$Combined.XCI.status[tab$Gene.name %in% merged.tab$geneSymbol])['escape']
  
 # inactive.biased<-inactive.biased/table(tab$Combined.XCI.status[tab$Gene.name %in% merged.tab$geneSymbol])['inactive']
  
  #variable.biased<-inactive.biased/table(tab$Combined.XCI.status[tab$Gene.name %in% merged.tab$geneSymbol])['variable']
  
  res.as<-rbind(res.as,c(paste(tissue.set,collapse='.'),escape.biased,inactive.biased,variable.biased))
}


colSums(apply(res.as[,2:4],2,as.numeric))

res.as<-rbind(res.as,c('Sum',colSums(apply(res.as[,2:4],2,as.numeric))))

#write.table(res.as,'/Users/karleg/Dimorph/XiASgenes.txt',sep='\t',row.names = F,quote = F)


##Significance

#DE
#wilcox.test(as.numeric(res.exp[1:(nrow(res.exp)-1),2]),as.numeric(res.exp[1:(nrow(res.exp)-1),3]),paired=TRUE)

#wilcox.test(as.numeric(res.exp[1:(nrow(res.exp)-1),2]),as.numeric(res.exp[1:(nrow(res.exp)-1),4]),paired=TRUE)


#AS
#wilcox.test(as.numeric(res.as[1:(nrow(res.as)-1),2]),as.numeric(res.as[1:(nrow(res.as)-1),3]),paired=TRUE)

#wilcox.test(as.numeric(res.as[1:(nrow(res.as)-1),2]),as.numeric(res.as[1:(nrow(res.as)-1),4]),paired=TRUE)



pval.exp<-fisher.test(rbind(as.integer(table(tab$Combined.XCI.status)),as.integer(res.exp[nrow(res.exp),2:4])))$p.value

pval.exp<-round(pval.exp,-log10(pval.exp)+3)

pval.as<-fisher.test(rbind(as.integer(table(tab$Combined.XCI.status)),as.integer(res.as[nrow(res.as),2:4])))$p.value

pval.as<-round(pval.as,-log10(pval.as)+3)



res.as<-as.data.frame(res.as)

res.as[,2:4]<-apply(res.as[,2:4],2,as.integer)

res.as[,2:4]<-res.as[,2:4]/rowSums(res.as[,2:4])

res.exp<-as.data.frame(res.exp)

res.exp[,2:4]<-apply(res.exp[,2:4],2,as.integer)

res.exp[,2:4]<-res.exp[,2:4]/rowSums(res.exp[,2:4])




df<-data.frame(Type=rep(sub('-.*','',colnames(res.as)[2:4]),3),data=c(rep('DE',3),rep('AS',3),rep('chrX',3)),count=as.numeric(c(res.exp[nrow(res.exp),2:4],res.as[nrow(res.as),2:4],table(tab$Combined.XCI.status)/sum(table(tab$Combined.XCI.status)))))


p<-ggplot(df, aes(data, count)) +  geom_bar(aes(fill = Type), position = "dodge", stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))


p<-p+ annotate("text", x=1, y=0.57, label=paste('P value:',pval.as ) )

p<-p+ annotate("text", x=3, y=0.72, label=paste('P value:',pval.exp ) )

p<-p+  ylab('Percent of Genes, %') +
  scale_y_continuous(labels = scales::percent)+xlab('')+ theme(text = element_text(size=14))+labs(fill="")

p

save.image('/Users/karleg/Dimorph/RDATA/figure2b.RData')


#Before running the following, use the Session menu to set working directory to source file location
load('figure2b.RData')
library(ggplot2)
library("ggsci")
library("gridExtra")


#NOTE pval.as=2.46e-20
#NOTE pval.exp=5.47e-50

p<-ggplot(df, aes(data, count)) +  geom_bar(aes(fill = Type), position = "dodge", stat="identity") + theme(axis.text.x = element_text(angle = 0, hjust = 1))
pn <- p + scale_fill_npg()  + theme_minimal() +  theme(text = element_text(size=20),axis.text.x = element_text(size=24, hjust=0.5),axis.title=element_blank())
pn<-pn+ annotate("text", x=1, y=0.57, label="italic(p)==2.46 %*% 10 ^ -20 ",parse=TRUE,size=6)
pn<-pn+ annotate("text", x=3, y=0.72, label="italic(p)==5.47 %*% 10 ^ -50 ", parse=TRUE,size=6)
pn +  ylab('') + scale_y_continuous(labels = scales::percent)+xlab('')+ theme(text = element_text(size=28))+labs(fill="")
pn

# cleaning the above code, and recommend figure style like this and saving figure this way to preserve text size
ggplot(df, aes(data, count)) +  
  geom_bar(aes(fill = Type), position = "dodge", stat="identity") + 
  annotate("text", x=1, y=0.57, label="italic(p)==2.46 %*% 10 ^ -20 ",parse=TRUE,size=3) +
  annotate("text", x=3, y=0.72, label="italic(p)==5.47 %*% 10 ^ -50 ", parse=TRUE,size=3) +
  ylab('') + xlab("")  +
  scale_fill_npg() + theme_bw() +
  theme(
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size=8, hjust=0.5),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 8)
        )
ggsave("figure2b_test.pdf", width = 4, height = 3)
 