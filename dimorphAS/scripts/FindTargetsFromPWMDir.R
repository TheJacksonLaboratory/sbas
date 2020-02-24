#This script finds RBP-target exon pair such that the RBP binds the target exon +-100bp,
#and models the psi values of the exon as a function of the expression levels of the RBP and sex

source("http://bioconductor.org/biocLite.R")

setwd('/projects/karleg/projects/dim')

library(biomaRt)

args <- commandArgs(trailingOnly = TRUE)

tissue<-'Adipose - Subcutaneous'

#This code runs on the cluster, each job handles 100 rMATs events.  The command line
#argument is the serial number of the first rMATs event to handle
rmats.file<-read.table('fromGTF.SE.txt',header=TRUE,skip=as.integer(args[1]),nrows=100) 

pwm.files<-list.files('RBP_PSSMs/')  #files with the PSSMs for the different RBP

#The table pwm.to.site contains the number of binding site for each RBP and rMATs event.
#The columns correspond to RBPs, and the rows to rMATs events.
pwm.to.site<-matrix(ncol=length(pwm.files),nrow=nrow(rmats.file)) 

#The following two function fill the table pwm.to.site with values
f1<-function(i)
{
  
  print(i)
  
  v=rmats.file[i,]  #For each row at index i of the rMATs events table
  
  lapply(1:length(pwm.files),f2,v,i) #For each RBP PSSM test if it binds the event v
}

f2<-function(j,v,i)
{
  c.num<-strsplit(as.character(v[[4]]),'chr')[[1]][2]
  
  #params will contain the path to the RBP PSSM file and the coordinates of the rMATS event
  
  params<-c(paste('RBP_PSSMs/',pwm.files[j],sep=''),paste(c.num,':',as.numeric(v[[6]]),'-',as.numeric(v[[7]]),sep=''),v[[1]])
  
  # The function run.tfbstools, which is implemented in TFBSToolsFromPWMDir.R, takes
  # params as argument and returns the set of binding sites of the RBP in the rMATs
  # exon 
  
  siteset<-run.tfbstools(params)
  
  #The function returns the number of binding sites
  length(siteset[pvalues(siteset, type="sampling")<=0.001]) 
}

#The following condition check if the binding matrix already exists, if it is not it
#creates it, otherwise is loads it from file.  This saves running time.
if (!file.exists(paste('binding_matrix',args[1],'.txt',sep='')))
{

	source('TFBSToolsFromPWMDir.R')

	res<-lapply(1:nrow(rmats.file),f1) #see documentation within f1 code above

	#The following call uses the function f1 to update the table pwm.to.site
	pwm.to.site<-do.call(rbind,lapply(res,unlist))

	#We next use biomaRt in order to map the RBP names from the binding matrices,
	#which are given as HGNC symbols, to ensembl IDs
	mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host="www.ensembl.org")

	hgnc.ensembl.map<-getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                        filters    = "hgnc_symbol",
                        values     = unlist(lapply(strsplit(pwm.files,'_'),'[',1)), 
                        mart       = mart)

	map.to.ens<-function(x,m){return (m[which(m[,2]==x)[1],1])}

	#Set the column names of pwm..to.site to the corresponding ensembl ids
	colnames(pwm.to.site)<-lapply(lapply(strsplit(pwm.files,'_'),'[',1),map.to.ens,hgnc.ensembl.map)

	pwm.to.site<-pwm.to.site[,!is.na(colnames(pwm.to.site))]

	rownames(pwm.to.site)<-rmats.file$ID

	write.table(pwm.to.site,paste('binding_matrix',args[1],'.txt',sep=''),sep='\t',quote = FALSE,row.names = TRUE,col.names = TRUE)
}else{
	
  #This is the case where the binding matrix already exists - in this case we simply laod it
	pwm.to.site<-read.table(paste('binding_matrix',args[1],'.txt',sep=''),sep='\t',row.names = 1,header = TRUE)

}

#load meta data for the GTEx dataset.
meta.data<-read.table('2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

#extract expression levels (counts) for all the RBPs.  We do that because the complete
#expressiob matrix is very large.
system(paste('cat GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct | grep -e Name -e ',paste(colnames(pwm.to.site),sep='',collapse = ' -e '),' > rbp_expression',args[1],'.txt',sep=''))

#restrict the meta data to the current tissue
meta.data<-meta.data[meta.data$body_site_s %in% tissue,]

#read the RBP expression levels that were extracted above
rbp.exp<-read.table(paste('rbp_expression',args[1],'.txt',sep=''),header=TRUE)

n.rbps<-nrow(rbp.exp) #number of RBP in the expressiobn data

#Next we format the tables rbp.exp and meta.data so that we can merge them, and will
#thus have a joint table of RBP expression and meta data 
rbp.exp<-t(rbp.exp)

colnames(rbp.exp)<-gsub("\\.[0-9]*","",rbp.exp[1,])

rbp.exp<-rbp.exp[-c(1,2),]

rownames(rbp.exp)=gsub('\\.','-',rownames(rbp.exp))

rownames(meta.data)<-meta.data$Sample_Name_s

rbp.exp<-merge(x=rbp.exp,y=meta.data,by='row.names')

#Now we need to get the psi values for the 100 events that this job handles:
psi.tab<-read.csv('rmats_final_matrix.se.jc.inc.txt',row.names = 1,header=TRUE,skip=as.integer(args[1]),nrows=100)

#The labels of the psi value correspond to the GTEx conditions/individuals.  We will use
#these labels to match the psi values of the rMATs events with the expression levels of
#the RBPs
psi.labels<-colnames(read.csv('rmats_final_matrix.se.jc.inc.txt',row.names = 1,header=TRUE,skip=0,nrows=1))

psi.labels<-gsub("\\.","-",psi.labels)

#Keep RBP expression data only for individuals who also have psi data
rbp.exp<-rbp.exp[rbp.exp$Row.names %in% psi.labels,]

#keep the labels of experiments who have psi data and RBP expression data
psi.labels.subset<-psi.labels[psi.labels %in% rbp.exp$Row.names]

#sanity check:
if (!(sum(rbp.exp$Row.names %in% psi.labels.subset)==nrow(rbp.exp) && sum( psi.labels.subset %in% rbp.exp$Row.names )==nrow(rbp.exp)))
  
  quit(status = 1)

#The following matrix will keep the p-values for the sex-RBP_expression interaction term 
#in the model
sex.mat<-matrix(nrow=nrow(rmats.file),ncol=n.rbps)

colnames(sex.mat)<-colnames(rbp.exp)[2:(1+n.rbps)]

rownames(sex.mat)<-rmats.file$ID

#The following matrix will keep the p-values for the RBP_expression term 
#in the model
rbp.mat<-matrix(nrow=nrow(rmats.file),ncol=n.rbps)

colnames(rbp.mat)<-colnames(rbp.exp)[2:(1+n.rbps)]

rownames(rbp.mat)<-rmats.file$ID

#The following two functions compute for each RBP and its targets the p-values for the
#RBP expression term and the sex-RBP_expression interaction term, but only for targets that
#are significantly different between males and females, as determined by the KS test
f4<-function(rbp,event.num,rbp.mat,sex.mat,psi.vals)
{
  
  if (pwm.to.site[event.num,rbp]==0)
    
    return()  #if thid RBP does not bind the target exon, go to the next RBP
  
  y=psi.vals #dependent variable
  
  #independent variable 1:  RBP expression
  x=as.numeric(as.character(rbp.exp[,which(colnames(rbp.exp)==rbp)]))
  
  #independent variable 2: RBP expression - sex interaction
  sex=ifelse(rbp.exp$sex_s=='male',1,0)
  
  sex.x<-sex*x
  
  df<-data.frame(x,y,sex.x)
  
  #fit quasi binomial model
  
  rbp.mod<-glm(y~x+sex.x,data=df,family=quasibinomial)
  
  s.mod<-summary(rbp.mod)
  
  #put results in the rbp and sex tables:
  
  rbp.mat[event.num,rbp]=s.mod$coefficients[2,4]
  
  sex.mat[event.num,rbp]=s.mod$coefficients[3,4]
  
  assign('rbp.mat',rbp.mat,envir=.GlobalEnv)
  
  assign('sex.mat',sex.mat,envir=.GlobalEnv)
  
}


f3<-function(event.num,psi.tab,rbp.mat,sex.mat)
{
  
  print(event.num) #The current rMATs event
  
  psi.vals<-as.numeric(psi.tab[event.num,]) #take the psi values for the current event
  
  psi.vals<-psi.vals[psi.labels %in% rbp.exp$Row.names] #Keep only the psi values that
  #have corresponding RBP expression values
  
  #order the psi values in the same order as the RBP expression values
  psi.vals<-psi.vals[order(match(psi.labels.subset,rbp.exp$Row.names))]

  #if the psi values if males and females are not significantly different, move to the next
  #rMATs event
  if (ks.test(psi.vals[rbp.exp$sex_s=='male'],psi.vals[rbp.exp$sex_s=='female'])$p.value>0.05)
    
    return()
  
  #If significantly different, check for each RBP that binds the exon (+-100bp) if there
  #is a significant correlation between RBP expression, sex and psi (done in the function f4)
  lapply(colnames(rbp.mat),f4,event.num,rbp.mat,sex.mat,psi.vals)
}

lapply(1:nrow(rmats.file),f3,psi.tab,rbp.mat,sex.mat)

write.table(sex.mat,paste('sex_interactions_p_value_matrix',args[1],'.txt',sep=''),sep='\t',quote = FALSE,row.names = TRUE,col.names = TRUE)

write.table(rbp.mat,paste('rbp_level_p_value_matrix',args[1],'.txt',sep=''),sep='\t',quote = FALSE,row.names = TRUE,col.names = TRUE)


