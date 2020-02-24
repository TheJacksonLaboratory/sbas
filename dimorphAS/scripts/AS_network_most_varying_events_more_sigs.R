library(coda)
library(rstan)
library(VGAM)

args <- commandArgs(trailingOnly = TRUE)

#Use these tissues because we need to have both sexually dimorphic events and non-sexually dimorphic events
#,and these are the tissues that have a relatively large number of sexually dimorphic events

tissue.list<-c('Heart - Left Ventricle',
               'Breast - Mammary Tissue',
               'Brain - Cortex.Brain - Frontal Cortex (BA9).Brain - Anterior cingulate cortex (BA24)',
               'Adrenal Gland'
               ,'Adipose - Subcutaneous','Muscle - Skeletal','Thyroid','Cells - Transformed fibroblasts','Artery - Aorta','Skin - Sun Exposed (Lower leg).Skin - Not Sun Exposed (Suprapubic)')


tissue=tissue.list[as.integer(args[1])]

events.table<-read.table(paste('other/',tissue,'se.txt',sep=''))

annot.table<-read.table('fromGTF.SE.txt',header=T)

rownames(annot.table)<-annot.table$ID

merged.table<-merge(events.table,annot.table,by='row.names')

all.genes<-read.table('GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct',sep='\t',header=T,skip=2,colClasses = c(rep("character", 2), rep("NULL", 11688)))

all.genes<-all.genes[!duplicated(all.genes$Description),]

#which(!merged.table$geneSymbol %in% all.genes$Description)

merged.table<-merged.table[order(merged.table$adj.P.Val),]

rbp.names<-unique(gsub('_.*','',list.files('RBP_PSSMs/')))

rbp.names<-rbp.names[rbp.names %in% all.genes$Description]

norm.counts<-read.table(paste('matrices/',tissue,'se_splice_levels.txt',sep=''))

norm.counts<-norm.counts[rownames(norm.counts) %in% merged.table$Row.names,]

norm.counts<-norm.counts[order(match(rownames(norm.counts),merged.table$Row.names)),]

norm.inc.counts<-norm.counts[,1:(ncol(norm.counts)/2)]

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

expression.mat<-read.table('GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct', nrow=1,sep='\t',header=T,skip=2)

colnames.expression.mat<-colnames(expression.mat)

meta.data<-read.table('2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

meta.data$Sample_Name_s<-gsub('-','\\.',meta.data$Sample_Name_s)

meta.data<-meta.data[meta.data$Sample_Name_s %in% colnames(expression.mat),]

col.in.tissue<-c()

for (col in colnames.expression.mat)
  
  col.in.tissue<-c(col.in.tissue, (col %in% meta.data$Sample_Name_s) && (meta.data$body_site_s[which(meta.data$Sample_Name_s==col)] %in% unlist(strsplit(tissue,split='\\.'))) && (meta.data$submitted_subject_id_s[which(meta.data$Sample_Name_s==col)]!='GTEX-11ILO'))


expression.mat<-read.table('GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct', colClasses = ifelse(col.in.tissue,"numeric","NULL"),sep='\t',header=T,skip=2)

sum(colnames(expression.mat) %in% colnames(norm.inc.counts))==ncol(expression.mat)

norm.inc.counts<-norm.inc.counts[,colnames(norm.inc.counts) %in% colnames(expression.mat)]

colnames(norm.skip.counts)<-gsub('\\.1$','',colnames(norm.skip.counts))

norm.skip.counts<-norm.skip.counts[,colnames(norm.skip.counts) %in% colnames(expression.mat)]

sum(colnames(norm.inc.counts) %in% colnames(expression.mat))==ncol(expression.mat)

sum(colnames(norm.skip.counts) %in% colnames(expression.mat))==ncol(expression.mat)

sum(colnames(norm.inc.counts)==colnames(norm.skip.counts))==ncol(norm.skip.counts)

expression.mat<-expression.mat[,order(match(colnames(expression.mat),colnames(norm.inc.counts)))]

meta.data<-meta.data[meta.data$Sample_Name_s %in% colnames(norm.inc.counts),]

meta.data<-meta.data[order(match(meta.data$Sample_Name_s,colnames(norm.inc.counts))),]

sum(colnames(expression.mat) == colnames(norm.inc.counts))==ncol(norm.inc.counts) #sanity check

all.genes<-read.table('GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct',sep='\t',header=T,skip=2,colClasses = c(rep("character", 2), rep("NULL", 11688)))

expression.mat<-expression.mat[!duplicated(all.genes$Description),]

all.genes<-all.genes[!duplicated(all.genes$Description),]

gene.names=unique(merged.table$geneSymbol)

which(!rbp.names %in% all.genes$Description)

expression.mat<-expression.mat[all.genes$Description %in% c(rbp.names,as.character(gene.names)),]

rownames.expression.mat<-all.genes$Description[all.genes$Description %in% c(rbp.names,as.character(gene.names))]

expression.mat<-expression.mat[!duplicated(rownames.expression.mat),]

rownames.expression.mat<-rownames.expression.mat[!duplicated(rownames.expression.mat)]

num.events=nrow(merged.table)

event.to.gene=c()

gexp=expression.mat[rownames.expression.mat %in% gene.names,]

rownames(gexp)<-rownames.expression.mat[rownames.expression.mat %in% gene.names]

gexp<-gexp[order(match(rownames(gexp),gene.names)),]

gexp=log2(gexp+0.5)

gexp=gexp-rowMeans(gexp)

gexp[apply(gexp,1,sd)>0,]=gexp[apply(gexp,1,sd)>0,]/apply(gexp[apply(gexp,1,sd)>0,],1,sd)

rexp=expression.mat[rownames.expression.mat %in% rbp.names,]

rownames(rexp)<-rownames.expression.mat[rownames.expression.mat %in% rbp.names]

rexp<-rexp[order(match(rownames(rexp),rbp.names)),]

rexp=log2(rexp+0.5)

rexp=rexp-rowMeans(rexp)

rexp=rexp/apply(rexp,1,sd)

for (i in (1:num.events))
  
  event.to.gene<-c(event.to.gene,which(unique(merged.table$geneSymbol)==merged.table[i,'geneSymbol']))

sex<-ifelse(meta.data$sex_s=='male',1,0)

sig<-ifelse(merged.table$adj.P.Val<=0.05 & abs(merged.table$logFC)>=log2(1.5),1,0)

write.table(merged.table,paste('EventsTable_',tissue,'.txt',sep='' ),sep='\t',quote = F,row.names=F)

dataList = list(
  as = round(norm.skip.counts) ,   #skip event counts across experiments
  c = round(norm.skip.counts+norm.inc.counts)    , #total counts for event, i.e. skip+inclusion, across experiments
  gexp = gexp, #read counts for genes (from gtex, take the raw counts) across experiments
  rexp = rexp, #read counts for RBPs (from gtex, take the raw counts)
  event_to_gene = event.to.gene,  #the gene index for each event (1 to the number of distinct genes) 
  Nrbp = length(rbp.names), #number of RBPs
  Nevents = nrow(merged.table),  #100 most varying AS events in 
  Nexp = ncol(expression.mat),#number of experiments such that we measured each event, gene and RBP in each experiment
  Ngenes = nrow(gexp),
  sex=sex,
  sig=sig
)


modelString = "
data {
int<lower=0> Nevents;
int<lower=0> Nexp;
int<lower=0> Nrbp;
int<lower=0> Ngenes;
int<lower=0> as[Nevents,Nexp] ;
int<lower=0> c[Nevents,Nexp] ;
matrix[Ngenes,Nexp] gexp ; 
matrix[Nrbp,Nexp] rexp ; 
int<lower=0> event_to_gene[Nevents];
int<lower=0,upper=1> sex[Nexp];
int<lower=0,upper=1> sig[Nevents];
}


parameters {
real beta0[Nevents] ;
real beta1[Nevents] ;
matrix[Nevents,Nrbp] beta2 ;
real beta3[Nevents];
matrix[Nrbp,2] beta4;

}
model {

for ( i in 1:Nexp ) {  


for ( j in 1:Nevents ) if (c[j,i]>0) { 

as[j,i] ~ binomial(c[j,i], inv_logit(beta0[j]+beta1[j]*sex[i]+dot_product(beta2[j,],rexp[,i])+beta3[j]*gexp[event_to_gene[j],i] ) );

}
}

for (k in 1:Nrbp){

for ( j in 1:Nevents ) { 

if (sig[j]==1)

  beta2[j,k] ~normal(beta4[k,1],1);

else

  beta2[j,k]~normal(beta4[k,2],1);

}
beta4[k,1]~normal(0,2);
beta4[k,2]~normal(0,2);
}


for ( j in 1:Nevents ) { 

beta1[j] ~ normal(0,2);
beta0[j] ~ normal(0,2);
beta3[j] ~ normal(0,2);
}

}
"

options(mc.cores = parallel::detectCores())

stanDso <- stan_model( model_code=modelString ) 

stanFit <- sampling( object=stanDso , data = dataList , chains = 3 ,iter = 5000,warmup=3000   , thin = 1,init=0 )

mcmcCoda = mcmc.list( lapply( 1:ncol(stanFit) , function(x) { mcmc(as.array(stanFit)[,x,]) } ) )

save( mcmcCoda , file=paste('McmcMostVaryingMoreSigs_',tissue,'.Rdata',sep='' ))




