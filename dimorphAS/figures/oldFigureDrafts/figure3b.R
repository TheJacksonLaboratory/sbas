
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

tissue <- tissue.list[[1]]
  
load(paste('/Users/karleg/Dimorph/McmcMostVaryingMoreSigs_',tissue,'.Rdata',sep=''))
  
mcmcCoda<-mcmcCoda[,which(grepl('beta2\\[101,87\\]',varnames(mcmcCoda))),drop=FALSE]

diagMCMC( mcmcCoda , parName=c("beta2[101,87]") )  

save.image('/Users/karleg/Dimorph/RDATA/figure3b.RData')

#Before running the following, use the Session menu to set working directory to source file location
load('figure3b.RData')
diagMCMC( mcmcCoda , parName=c("beta2[101,87]") )  

