library(limma)

library(multtest)

library(Biobase)

library(edgeR)

#DS analysis

inc.iso.counts.mem<-read.csv('/Users/karleg/Downloads/rmats_final.se.jc.ijc.txt',header=TRUE)

skip.iso.counts.mem<-read.csv('/Users/karleg/Downloads/rmats_final.se.jc.sjc.txt',header=TRUE)

meta.data<-read.csv('/Users/karleg/Downloads/SraRunTable.noCram.noExome.noWGS.totalRNA.txt',header=TRUE)

meta.data<-meta.data[meta.data$Run %in% colnames(inc.iso.counts.mem),]

rownames(inc.iso.counts.mem)<-inc.iso.counts.mem$ID

inc.iso.counts.mem<-inc.iso.counts.mem[,-1]

inc.iso.counts.mem<-inc.iso.counts.mem[,order(match(colnames(inc.iso.counts.mem),meta.data$Run))]

inc.iso.counts.mem<-inc.iso.counts.mem[,!grepl('11ILO',meta.data$Sample.Name) & (meta.data$body_site %in% 'Breast - Mammary Tissue')]

rownames(skip.iso.counts.mem)<-skip.iso.counts.mem$ID

skip.iso.counts.mem<-skip.iso.counts.mem[,-1]

skip.iso.counts.mem<-skip.iso.counts.mem[,order(match(colnames(skip.iso.counts.mem),meta.data$Run))]

skip.iso.counts.mem<-skip.iso.counts.mem[,!grepl('11ILO',meta.data$Sample.Name) & (meta.data$body_site %in% 'Breast - Mammary Tissue')]

meta.data<-meta.data[!grepl('11ILO',meta.data$Sample.Name),]

meta.data<-meta.data[meta.data$body_site %in% 'Breast - Mammary Tissue',]

keep.cols<-colSums(inc.iso.counts.mem>0)>=0.25*nrow(inc.iso.counts.mem)

keep.cols<-keep.cols & colSums(skip.iso.counts.mem>0)>=0.25*nrow(skip.iso.counts.mem)

inc.iso.counts.mem=inc.iso.counts.mem[,keep.cols]

skip.iso.counts.mem=skip.iso.counts.mem[,keep.cols]

meta.data=meta.data[keep.cols,]

if (sum(colnames(skip.iso.counts.mem)==meta.data$Run)!=nrow(meta.data))
  
  quit(1)

if (sum(colnames(inc.iso.counts.mem)==meta.data$Run)!=nrow(meta.data))
  
  quit(1)

exprs.mat <- matrix(data=cbind(as.matrix(skip.iso.counts.mem),as.matrix(inc.iso.counts.mem)),ncol=2*ncol(skip.iso.counts.mem),nrow=nrow(skip.iso.counts.mem),dimnames=list(1:nrow(skip.iso.counts.mem),c(colnames(skip.iso.counts.mem),paste(colnames(skip.iso.counts.mem),'2',sep=''))))

keep.events<-rep(T,nrow(exprs.mat))

groups=c(paste0(meta.data$sex,'-skip'),paste0(meta.data$sex,'-inc'))

groups<-groups[colSums(exprs.mat)>0]

exprs.mat<-exprs.mat[,colSums(exprs.mat)>0]

zero.rows=rep(F,nrow(exprs.mat))

for (group in c('male-skip','female-skip','male-inc','female-inc'))
  
  zero.rows=zero.rows | (rowSums(exprs.mat[,groups %in% group] >0) < (1/4)*min(table(groups)))

exprs.mat[zero.rows,]=0

for (group in c('male-skip','female-skip','male-inc','female-inc'))
  
  keep.events<-keep.events & (rowSums((exprs.mat[,groups %in% group]) >= 1) >= (1/4)*min(table(groups)))

rm(exprs.mat)

inc.counts.male<-inc.iso.counts.mem[,meta.data$sex=='male']

inc.counts.female<-inc.iso.counts.mem[,meta.data$sex=='female']

skip.counts.male<-skip.iso.counts.mem[,meta.data$sex=='male']

skip.counts.female<-skip.iso.counts.mem[,meta.data$sex=='female']

inc.counts.mat<-cbind(inc.counts.male,inc.counts.female)

rownames(inc.counts.mat)<-rownames(inc.iso.counts.mem)

skip.counts.mat<-cbind(skip.counts.male,skip.counts.female)

rownames(skip.counts.mat)<-rownames(skip.iso.counts.mem)

counts.mat<-cbind(inc.counts.male,inc.counts.female,skip.counts.male,skip.counts.female)

rownames(counts.mat)<-rownames(inc.iso.counts.mem)

counts.mat<-counts.mat[keep.events,]

isoform<-c(rep(1,ncol(inc.counts.male)+ncol(inc.counts.female)),rep(0,ncol(skip.counts.male)+ncol(skip.counts.female)))

sex<-c(rep(1,ncol(inc.counts.male)),rep(0,ncol(inc.counts.female)),rep(1,ncol(skip.counts.male)),rep(0,ncol(skip.counts.female)))

block=rep(c(1:ncol(inc.counts.male),(1+ncol(inc.counts.male)):(ncol(inc.counts.female)+ncol(inc.counts.male))),2)

keep.cols<-colSums(counts.mat)>0

isoform<-isoform[keep.cols]

sex<-sex[keep.cols]

block<-block[keep.cols]

counts.mat<-counts.mat[,keep.cols]

norm.factors=calcNormFactors(DGEList(counts=counts.mat[,isoform==1]+counts.mat[,isoform==0]))

counts.mat=DGEList(counts=counts.mat)

counts.mat$samples$norm.factors=rep(norm.factors$samples$norm.factors,2)

design <- model.matrix(~sex+isoform+sex*isoform)

correlation =  duplicateCorrelation(voom(counts.mat)$E,design,block=block)$consensus.correlation

v <- voom(counts.mat, design,block = block, correlation = correlation)

fit <- lmFit(v, design,block = block, correlation = correlation)

fit <- eBayes(fit, robust=TRUE)

res=topTable(fit, coef='sex:isoform',number=nrow(counts.mat))

write.table(res[res$adj.P.Val<=0.05 & abs(res$logFC)>=log2(1.5),],paste(paste('Breast - Mammary Tissue',collapse='.'),'se.txt',sep=''),row.names = T,col.names = T,quote = F)

write.table(rownames(res),paste0(paste('Breast - Mammary Tissue',collapse='.'),'all_genes.txt'))


#Run Ontologizer

meta.data<-read.table('/Users/karleg/Downloads/fromGTF.SE.txt',sep='\t',header=TRUE)

de.tab<-read.table('Breast - Mammary Tissuese.txt')

de.tab.with.meta<-merge(de.tab,meta.data,by.x='row.names',by.y='ID')

all.genes<-read.table('Breast - Mammary Tissueall_genes.txt')

write.table(meta.data$geneSymbol[meta.data$ID %in% all.genes[,1]],'universe.txt',quote = F,row.names = F,col.names = F)

write.table(de.tab.with.meta$geneSymbol,'gene_set.txt',quote = F,row.names = F,col.names = F)

system('java -jar /Users/karleg/Ontologizer/Ontologizer.jar -g /Users/karleg/Ontologizer/go.obo -a /Users/karleg/Ontologizer/goa_human.gaf -s gene_set.txt -p universe.txt -c Term-For-Term -m Benjamini-Hochberg -n')

