
event.type='SE'

univ.tab<-NULL

max.rows<-(round(nrow(read.table(paste0('/Users/karleg/Dimorph/fromGTF.',event.type,'.txt')))/1000))*1000

for (i in seq(0,max.rows,1000))
  
  univ.tab<-rbind(univ.tab,read.csv(paste('/Users/karleg/Dimorph/domain_plots_all/',event.type,'found_annotation_all_genes',i,'.txt',sep=''),sep='\t'))

write.table(univ.tab,paste0('/Users/karleg/Dimorph/domain_plots_all/',event.type,'found_annotation_all_genes.txt'),quote = F,row.names = F,col.names = T,sep='\t')



setwd('/Users/karleg/Ontologizer/')

meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

tissue.sets<-c('Heart - Left Ventricle',
               'Breast - Mammary Tissue',
               'Brain - Cortex.Brain - Frontal Cortex (BA9).Brain - Anterior cingulate cortex (BA24)',
               'Adrenal Gland'
               ,'Adipose - Subcutaneous','Muscle - Skeletal','Thyroid','Cells - Transformed fibroblasts','Artery - Aorta','Skin - Sun Exposed (Lower leg).Skin - Not Sun Exposed (Suprapubic)')

is.dimorph='Yes'

affected.by.expression=TRUE

test.mat<-matrix(ncol=45,nrow=2)

test.mat[,1]<-c('Peptidase','domain')

test.mat[,2]<-c('Proteinkinase','domain')

test.mat[,3]<-c('Acyl-CoA','domain')

test.mat[,4]<-c('phosphothreonine','mod.res')

test.mat[,5]<-c('Thioredoxin','domain')

test.mat[,6]<-c('NAD-binding','domain')

test.mat[,7]<-c('phosphoserine','mod.res')

test.mat[,8]<-c('helical','transmem')

test.mat[,9]<-c('TPR_REGION','domain')

test.mat[,10]<-c('Phosphotyrosine','mod.res')

test.mat[,11]<-c('Nuclearlocalizationsignal','motif')

test.mat[,12]<-c('N-acetylmethionine','mod.res')

test.mat[,13]<-c('EF-hand','domain')

test.mat[,14]<-c('RING-type','domain')

test.mat[,15]<-c('Ubiquitin-like','domain')

test.mat[,16]<-c('Homeobox','domain')

test.mat[,17]<-c('Protonacceptor','act-site')

test.mat[,18]<-c('Cytoplasmic','topo-dom')

test.mat[,19]<-c('Lumenal','topo-dom')

test.mat[,20]<-c('Mitochondrialmatrix','topo-dom')

test.mat[,21]<-c('Mitochondrialinnermembrane','intramem')

test.mat[,22]<-c('N6-acetyllysine','mod.res')

test.mat[,23]<-c('Dynamin-typeG','domain')

test.mat[,24]<-c('Extracellular','topo-dom')

test.mat[,25]<-c('Adaptin','domain')

test.mat[,26]<-c('Bbox','domain')

test.mat[,27]<-c('Hist_deacetyl','domain')

test.mat[,28]<-c('WD_REPEATS_REGION','domain')

test.mat[,29]<-c('BTB','domain')

test.mat[,30]<-c('EGF-like','domain')

test.mat[,31]<-c('Fibronectintype','domain')

test.mat[,32]<-c('Calponin-homology','domain')

test.mat[,33]<-c('Ion_trans','domain')

test.mat[,34]<-c('Ig-like','domain')

test.mat[,35]<-c('Lyase','domain')

test.mat[,36]<-c('ANK_REP_REGION','domain')

test.mat[,37]<-c('phosphatase','domain')

test.mat[,38]<-c('Myosinmotor','domain')

test.mat[,39]<-c('Ras-GEF','domain')

test.mat[,40]<-c('ATP-Binding','domain')

test.mat[,41]<-c('histone','domain')

test.mat[,42]<-c('methyltransferase','domain')

test.mat[,43]<-c('tubulin','domain')

test.mat[,44]<-c('NAD-binding_dom','domain')

test.mat[,45]<-c('Myosinmotor','domain')

res<-matrix(nrow=ncol(test.mat),ncol=6)

colnames(res)<-c('p-value','Domain','UniverseTotal','UniverseDomain','NMDTotal','NMDDomain')

for (i in (1:ncol(test.mat)))
{
  
  mod.type=test.mat[1,i]
  
  domain.field=test.mat[2,i]  #act.site , domain, motif , binding, topo.dom  etc.  - colnames(domain.tab) gives the full list
  
  nmd.tab<-read.table('/Users/karleg/Dimorph/NMD_summary.txt',sep='\t')
  
  universe.tab<-read.csv(paste0('/Users/karleg/Dimorph/domain_plots_all/',event.type,'found_annotation_all_genes.txt'),header=T,sep='\t')
  
  universe.tab<-universe.tab[universe.tab$ID %in% (rownames(nmd.tab)[as.integer(nmd.tab[,2])>0]),]
  
  dom.field.num<-which(colnames(universe.tab)==domain.field)

  univ.total<-length(unique(universe.tab$ID)) 
  
  univ.domain<-length(unique(universe.tab[grepl(tolower(mod.type),tolower(as.character(universe.tab[,dom.field.num]))),'ID']))
   
  nmd.tab<-nmd.tab[as.integer(nmd.tab[,2])>0,]
  
  nmd.total<-sum(nmd.tab$num.nmd>0)
  
  nmd.domain<-universe.tab[grepl(tolower(mod.type),tolower(as.character(universe.tab[,dom.field.num]))),]
  
  nmd.domain<-nmd.domain[nmd.domain$ID %in% rownames(nmd.tab[nmd.tab[,1]>0,]),]
  
  nmd.domain<-length(unique(nmd.domain$ID))
  
  res[i,1]<-phyper(q=nmd.domain, m=univ.domain, n=univ.total-univ.domain, k=nmd.total)
  
  res[i,2]<-mod.type
  
  res[i,3]<-univ.total
  
  res[i,4]<-univ.domain
  
  res[i,5]<-nmd.total
  
  res[i,6]<-nmd.domain
  
}

res<-cbind(res,rep(0,nrow(res)))

colnames(res)[ncol(res)]<-'FDR'

res[,ncol(res)]<-mt.rawp2adjp(as.numeric(res[,1]),proc='BH')[[1]][,2][order(mt.rawp2adjp(as.numeric(res[,1]),proc='BH')$index)]

res<-res[order(as.numeric(res[,1])),]

write.table(res[,c(2,3,4,5,6,1,7)],'NMD_Domain_Depletion.txt',sep='\t',quote = F,row.names = F)

plot.tab<-read.table('/Users/karleg/Ontologizer/NMD_Domain_Depletion.txt',sep='\t',header = T)

plot.tab$Domain<-as.character(plot.tab$Domain)

plot.tab$UniverseTotal<-plot.t

df<-data.frame(Type=rep(colnames(plot.tab)[2:5],sum(plot.tab$FDR<=0.05)),Domain=rep(unlist(lapply(plot.tab$Domain[plot.tab$FDR<=0.05],rep,4)),sum(plot.tab$FDR<=0.05)),
               Count=(as.integer(unlist(split(plot.tab[plot.tab$FDR<=0.05,2:5],1:sum(plot.tab$FDR<=0.05))))))

df<-df[!duplicated(paste(df$Type,df$Domain)),]

df2<-matrix(ncol=3,nrow=0)

colnames(df2)<-colnames(df)

for (i in seq(1,nrow(df),4))
{
  
  df2<-rbind(df2,c('%Universe',as.character(df$Domain)[i],df$Count[i+1]/df$Count[i]))
  
  df2<-rbind(df2,c('%NMD',as.character(df$Domain)[i],df$Count[i+3]/df$Count[i+2]))
  
}

df2<-as.data.frame(df2)

df2$Count<-as.numeric(df2$Count)

library(ggplot2)
library(ggsci)

p<-ggplot(df2, aes(Domain, Count)) +  geom_bar(aes(fill = Type), position = "dodge", stat="identity")+ylab('Percentage') + theme(axis.text.x = element_text(angle = 90, hjust = 1),text = element_text(size=20))+
  scale_fill_npg()  

for (i in 1:sum(plot.tab$FDR<=0.05))

  p<-p+ annotate("text", x=i, y=15.5, label=paste('P:',round(plot.tab$FDR[i],-log10(plot.tab$FDR[i])+3) ) )

p                 

save.image('/Users/karleg/Dimorph/RDATA/figure4e.RData')   

load('figure4e.RData')


library(ggplot2)
library(ggsci)


p<-ggplot(df2, aes(Domain, Count)) +  geom_bar(aes(fill = Type), position = "dodge", stat="identity")+ylab('Percentage') + theme_minimal()
p<- p+theme(axis.text.x = element_text(size=24,angle = 90, hjust = 1),
            text = element_text(size=20),
            axis.text.y = element_text(size=18),
            axis.title.x = element_blank(), 
            legend.title = element_blank())
p <- p+ scale_fill_npg()  

for (i in 1:sum(plot.tab$FDR<=0.05)) {
  if (plot.tab$FDR[i]<10^-20) {
    p<-p+ annotate("text", x=i, y=15.5, label="***",size=8) 
  } else if (plot.tab$FDR[i]<10^-5) {
    p<-p+ annotate("text", x=i, y=15.5, label="**",size=8) 
  } else if (plot.tab$FDR[i]<0.05) {
    p<-p+ annotate("text", x=i, y=15.5, label="*",size=8) 
  } 
  #p<-p+ annotate("text", x=i, y=15.5, label=paste('P:',round(plot.tab$FDR[i],-log10(plot.tab$FDR[i])+3) ) )
}
p                 

                     
                     
   