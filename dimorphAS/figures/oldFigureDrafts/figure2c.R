


print.pathway<-function(tissue.set)
{
  
  gene.associations<-read.table('anno-gene_set-MGSA-Benjamini-Hochberg.txt',sep='\t')
  
  sig.categories<-read.table('table-gene_set-MGSA-Benjamini-Hochberg.txt',sep='\t',header=T)
  
  sig.categories<-sig.categories[sig.categories$marg>=0.5,c('ID','name')]
  
  if (nrow(sig.categories)>0)
    
    for (i in (1:nrow(sig.categories)))
    {
      
      cat.genes<-as.character(gene.associations$V1[grepl(sig.categories[i,1],gene.associations$V3)])
      
      write.table(cat.genes,paste('/Users/karleg/Ontologizer/pathway_genes/',paste(tissue.set,collapse='.'),'.',gsub('/','',sig.categories[i,2]),sep=''),row.names = F,col.names = F,quote = F)
      
    }
  
  
}



library(multtest)

source('/Users/karleg/Dimorph/combine_hits.R')

setwd('/Users/karleg/Ontologizer/')

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


joint.tab<-matrix(nrow=0,ncol=5)

colnames(joint.tab)<-c('tissue','marg','Term','Event','Test')

for (ignore.nmd in c(T,F))
  for (run.on.nmd in c(T,F))
    for (tissue.set in tissue.sets)
    {
      for (event.type in c('se','mxe','ri','a3ss','a5ss')){
        
        if (ignore.nmd==F && event.type!='se')
          
          next
  
       if (!file.exists(paste0('/Users/karleg/Dimorph/other/',paste(tissue.set,collapse='.'),event.type,'.txt')))
    
        next
  
      de.tab<-read.table(paste0('/Users/karleg/Dimorph/other/',paste(tissue.set,collapse='.'),event.type,'.txt'))
      
      an.tab<-read.table(paste0('/Users/karleg/Dimorph/fromGTF.',toupper(event.type),'.txt'),header=T)
      
      nmd.tab<-read.table('/Users/karleg/Dimorph/NMD_summary.txt',sep='\t')
    
      test.string='Ignore NMD'
        
      if (event.type=='se' && !ignore.nmd && run.on.nmd)
      {
        an.tab<-an.tab[as.integer(nmd.tab[,2])>0 & nmd.tab[,1]==0,]
        
        test.string='Non-NMD Events'
      
      }else if (event.type=='se' && !ignore.nmd ){
        
        an.tab<-an.tab[as.integer(nmd.tab[,2])>0 & nmd.tab[,1]>0,]
        
        test.string='NMD Events'
        
      }
      
      merged.tab<-merge(an.tab,de.tab,by='row.names')
      
      de.genes<-as.character(unique(merged.tab$geneSymbol[abs(merged.tab$logFC)>=log2(1.5) & merged.tab$adj.P.Val<=0.05]))
      
      all.genes<-as.character(unique(merged.tab$geneSymbol))
      
      
      if (length(de.genes)<5)
        
        next
      
      write.table(de.genes,'/Users/karleg/Ontologizer/gene_set.txt',quote = F,row.names = F,col.names = F)
      
      write.table(all.genes,'/Users/karleg/Ontologizer/universe_e.txt',quote = F,row.names = F,col.names = F)
      
      system('java -jar Ontologizer.jar -g go.obo -a goa_human.gaf -s gene_set.txt -p universe_e.txt -c MGSA -m Benjamini-Hochberg -n')
      
      ont.out<-read.table('table-gene_set-MGSA-Benjamini-Hochberg.txt',sep='\t',header=T)
      
    #  print.pathway(tissue.set)
      
      next.tab<-matrix(c(rep(paste(tissue.set,collapse='.'),length(ont.out$marg)),ont.out$marg,as.character(ont.out$name),rep(event.type,length(ont.out$marg)),rep(test.string,length(ont.out$marg))),ncol=5)
      
      joint.tab<-rbind(joint.tab,next.tab[as.numeric(next.tab[,2])>0.5,])
      
  }
}


View(table(joint.tab[,1]))

View(table(joint.tab[,3]))

table(table(joint.tab[,3]))

which(table(joint.tab[,3])==max(table(joint.tab[,3])))

write.table(joint.tab,'GO-enrichment_AS_MGSA.txt',sep='\t',quote = F)

rows.to.plot<-c()

rows.to.plot<-c(rows.to.plot,which(joint.tab[,1]=='Heart - Left Ventricle' & joint.tab[,3]=='alpha-actinin binding'))

#rows.to.plot<-c(rows.to.plot,which(joint.tab[,1]=='Heart - Left Ventricle' & joint.tab[,3]=='cellular lipid metabolic process'))

#rows.to.plot<-c(rows.to.plot,which(joint.tab[,1]=='Breast - Mammary Tissue'& joint.tab[,3]=='proteinaceous extracellular matrix'))

#rows.to.plot<-c(rows.to.plot,which(joint.tab[,1]=='Breast - Mammary Tissue'& joint.tab[,3]=='semaphorin receptor binding'))

#get AS fold-changes

#get gene expression fold-changes

pway.as<-c()

data.de<-c()

event.as<-c()

pway.de<-c()

data.as<-c()

tissue.de<-c()



for (i in rows.to.plot)
{
  next.genes<-as.character(read.table(paste('/Users/karleg/Ontologizer/pathway_genes/',joint.tab[i,1],'.',gsub('/','',joint.tab[i,3]),sep=''))[,1])
  
  fold.changes.as<-read.table(paste0('/Users/karleg/Dimorph/other/',joint.tab[i,1],joint.tab[i,4],'.txt'))
  
  an.tab<-read.table(paste0('/Users/karleg/Dimorph/fromGTF.',joint.tab[i,4],'.txt'),header=T)
  
  merged.tab<-merge(an.tab,fold.changes.as,by='row.names')
  
  data.as<-c(data.as,merged.tab$logFC[(merged.tab$geneSymbol %in% next.genes) & (merged.tab$adj.P.Val<=0.05) & (abs(merged.tab$logFC)>=log2(1.5))])
  
  event.as<-c(event.as,rep(toupper(joint.tab[i,4]),sum((merged.tab$geneSymbol %in% next.genes) & (merged.tab$adj.P.Val<=0.05) & (abs(merged.tab$logFC)>=log2(1.5)))))
  
  pway.as<-c(pway.as,rep(joint.tab[i,3],sum((merged.tab$geneSymbol %in% next.genes) & (merged.tab$adj.P.Val<=0.05) & (abs(merged.tab$logFC)>=log2(1.5)))))
    
  fold.changes.de<-read.table(paste0('/Users/karleg/Dimorph/gene_expression/all_genes/DE_result_',joint.tab[i,1],'.txt'))
  
  fold.changes.de<-fold.changes.de[rownames(fold.changes.de) %in% next.genes,]
  
  pway.de<-c(pway.de,rep(joint.tab[i,3],sum(rownames(fold.changes.de) %in% next.genes)))
  
  data.de<-c(data.de,fold.changes.de$logFC)
  
  tissue.de<-c(tissue.de,rep(joint.tab[i,1],sum(rownames(fold.changes.de) %in% next.genes)))
  
}

library(ggplot2)

pway <- unique(pway.de)[[1]]

df<-data.frame(logFC=c(data.de[pway.de==pway],data.as[pway.as==pway]),type=c(rep('Expression',sum(pway.de==pway)),rep("Inclusion",sum(pway.as==pway))))

p <-  ggplot(df, aes(x=type, y=logFC,color=type)) + 
  geom_boxplot()+ggtitle(paste0(pway,', ',unique(tissue.de[pway.de==pway]),'(',unique(event.as[pway.as==pway]),')'))+
  ylab('logFC')+xlab('')+ theme(axis.ticks.y=element_blank(),axis.text.y=element_blank())+
  coord_flip()+geom_hline(yintercept=0,color='red',linetype='dashed')+
  theme(axis.text=element_text(size=16))+
  theme(axis.title=element_text(face="bold",size="14",color="black"),axis.text=element_text(size=14,face="bold"))+
  theme(plot.title = element_text(color="black", size=14, face="bold")) 

legend = gtable_filter(ggplotGrob(p), "guide-box") 

p<-list()

for (pway in unique(pway.de))
{

df<-data.frame(logFC=c(data.de[pway.de==pway],data.as[pway.as==pway]),type=c(rep('Expression',sum(pway.de==pway)),rep("Inclusion",sum(pway.as==pway))))
  
p[[length(p)+1]] <- ggplot(df, aes(x=type, y=logFC,color=type)) + geom_boxplot()+ggtitle(paste0(pway,', ',unique(tissue.de[pway.de==pway]),'(',unique(event.as[pway.as==pway]),')'))+ylab('')+xlab('')+ theme(legend.position="none",axis.ticks.y=element_blank(),axis.text.y=element_blank())+ coord_flip()+geom_hline(yintercept=0,color='red',linetype='dashed')+theme(axis.text=element_text(size=16))+theme(plot.background = element_rect(size=0.5,linetype="solid",color="black"))


}

p[[1]]

save.image('/Users/karleg/Dimorph/RDATA/figure2c.RData')

#grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]], nrow = 3,
      #   3    bottom = textGrob(bquote('logFC'),vjust=-25, gp = gpar(fontface = "bold", cex = 1.5)),legend)

#Before running the following, use the Session menu to set working directory to source file location
load('figure2c.RData')
library(ggplot2)
pn<-ggplot(df,aes(factor(type),logFC)) + geom_violin(aes(fill=type)) + scale_fill_npg() 
pn <- pn + theme_minimal() +  theme(text = element_text(size=20),
                                    axis.text = element_text(size=24, hjust=0.5),
                                    axis.title.x=element_blank(),
                                    axis.title.y = element_text(size=24),
                                    plot.title = element_text(hjust = 0.5),
                                    legend.position = "none") + ylab("log fold change") + labs(title="alpha-actinin binding")
pn <- pn +  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
pn




