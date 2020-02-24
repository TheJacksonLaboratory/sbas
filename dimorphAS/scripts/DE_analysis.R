setwd('/Users/karleg/Dimorph/')

meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

our.genes<-c()

for (event.type in c('se','mxe','ri','a3ss','a5ss'))
{
  
  events.tab<-read.table(paste('/Users/karleg/Dimorph/fromGTF.',toupper(event.type),'.txt',sep=''),header=T)
  
  our.genes<-unique(c(our.genes,as.character(events.tab$geneSymbol)))
  
}


tissue.sets<-unique(as.character(meta.data$body_site_s))
  
for (tissue.set in tissue.sets)
{
  
expression.mat<-read.table('GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct', nrow=1,sep='\t',header=T,skip=2)
  
colnames.expression.mat<-colnames(expression.mat)

meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

if (length(table(meta.data$sex_s[meta.data$body_site_s %in% tissue.set]))==1)
  
  next

meta.data$Sample_Name_s<-gsub('-','\\.',meta.data$Sample_Name_s)

meta.data<-meta.data[meta.data$Sample_Name_s %in% colnames(expression.mat),]

col.in.tissue<-c()

if (sum(meta.data$body_site_s %in% tissue.set)==0)
  
  next

for (col in colnames.expression.mat)

  col.in.tissue<-c(col.in.tissue, (col %in% meta.data$Sample_Name_s) && (meta.data$body_site_s[which(meta.data$Sample_Name_s==col)] %in% tissue.set))

expression.mat<-read.table('GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct', colClasses = ifelse(col.in.tissue,"numeric","NULL"),sep='\t',header=T,skip=2)
  
write.table(expression.mat,paste('gene_expression/expression',tissue.set,'.txt',sep=''),sep='\t',quote = F)

}


for (tissue.set in tissue.sets)
{

all.genes<-read.table('GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct',sep='\t',header=T,skip=2,colClasses = c(rep("character", 2), rep("NULL", 11688)))

if (file.size(paste('gene_expression/expression',tissue.set,'.txt',sep=''))<=1)
  
  next

x <- read.delim(paste('gene_expression/expression',tissue.set,'.txt',sep=''))

#x<-x[all.genes$Description %in% our.genes,]

#all.genes<-all.genes[all.genes$Description %in% our.genes,2]

x<-x[!duplicated(all.genes[,'Description']),]

all.genes<-all.genes[!duplicated(all.genes[,'Description']),'Description']
  
rownames(x)<-all.genes

#meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

#meta.data$Sample_Name_s<-gsub('-','\\.',meta.data$Sample_Name_s)

#meta.data<-meta.data[meta.data$Sample_Name_s %in% colnames(x),]

pheno<-read.csv('GTEx_v7_Annotations_SubjectPhenotypesDS.txt',sep='\t')

pheno$SUBJID<-gsub('-','\\.',pheno$SUBJID)

sex<-c()  #the value 2 in the phenotypic data is the values that is 1 in design.  This value correslonds to female.

age<-c()

cod<-c()

for (col in colnames(x))
{
  sex<-c(sex,pheno$SEX[which(pheno$SUBJID==paste(unlist(strsplit(col,'\\.'))[c(1,2)],collapse='.'))])
  
  age<-c(age,pheno$AGE[which(pheno$SUBJID==paste(unlist(strsplit(col,'\\.'))[c(1,2)],collapse='.'))])
  
  cod<-c(cod,pheno$DTHHRDY[which(pheno$SUBJID==paste(unlist(strsplit(col,'\\.'))[c(1,2)],collapse='.'))])
  
}

if (length(table(sex))==1)
  
  next

cod[is.na(cod)]=0

#group <- meta.data$sex_s

y <- DGEList(counts=x,group=factor(sex))

y <- calcNormFactors(y)

groups<-sex

keep.events<-rep(T,nrow(y))

for (group in c(1,2))
  
  keep.events<-keep.events & (rowSums(cpm(y[,groups %in% group]) > 1) >= 0.25*min(table(groups)))

if (sum(keep.events)==0)
  
  next

y<-y[keep.events,]  #if only using the logFC to compare with AS then do not screen

design <- model.matrix(~factor(sex))

v <- voom(y, design)

fit <- lmFit(v, design)

fit <- eBayes(fit, robust=TRUE)

res=topTable(fit, coef='factor(sex)2',number=nrow(y))

write.table(res,paste('gene_expression/all_genes/DE_result_',tissue.set,'.txt',sep=''),sep='\t',quote = F)

}
#edgeR

#y <- estimateDisp(y,design)

#fit <- glmFit(y,design)

#lrt <- glmLRT(fit,coef='factor(sex)2')

#hist(lrt$table$PValue)




