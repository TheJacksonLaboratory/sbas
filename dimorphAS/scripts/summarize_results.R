#This script assembles the result matrices of FindTargetsFromPWMDir.R that were created on the
#cluster, and plots the results

tissue<-'Adipose - Subcutaneous'

level.tab<-matrix(ncol=86,nrow=0)

#read in the RBP expression term p-values matrix:
for (i in seq(0,39100,100))
{
  if (file.exists(paste('/Users/karleg/Dimorph/',tissue,'/rbp/rbp_level_p_value_matrix',i,'.txt',sep='')))
  {  
    level.tab<-rbind(level.tab,read.table(paste('/Users/karleg/Dimorph/',tissue,'/rbp/rbp_level_p_value_matrix',i,'.txt',sep=''),sep='\t'))
  }else{
    print(paste('rbp_level_p_value_matrix',i,'.txt does not exist',sep=''))
    
    level.tab<-rbind(level.tab,matrix(rep(NA,86*100),ncol=86,dimnames=list((i+1):(i+100),colnames(level.tab))))
  }  
}

#read in the RBP_expression sex interaction term p-values matrix:
sex.tab<-matrix(ncol=86,nrow=0)

for (i in seq(0,39100,100))
{
  if (file.exists(paste('/Users/karleg/Dimorph/',tissue,'/sex/sex_interactions_p_value_matrix',i,'.txt',sep='')))
  { 
    sex.tab<-rbind(sex.tab,read.table(paste('/Users/karleg/Dimorph/',tissue,'/sex/sex_interactions_p_value_matrix',i,'.txt',sep=''),sep='\t'))
  }else{
    print(paste('sex_interactions_p_value_matrix',i,'.txt does not exist',sep=''))
    
    sex.tab<-rbind(sex.tab,matrix(rep(NA,86*100),ncol=86,dimnames=list((i+1):(i+100),colnames(level.tab))))
    
  } 
}

#For each RBP and its target such that rthe target has differential psi between males and
#females, collect the p-values of the model into the matrix 'results'
results=matrix(nrow=0,ncol=2)

colnames(results)<-c('level','sex')

for (i in (1:39180))
{
  
  for (j in (1:86))
  {
  
    if (!is.na(level.tab[i,j] | !is.na(sex.tab[i,j])))
    {
      
      results<-rbind(results,c(level.tab[i,j],sex.tab[i,j]))
      
      rownames(results)[nrow(results)]=paste(colnames(level.tab)[j],i,sep='-')
      
    }
      
  }
  
}


#This function creates a scatter plot of RBP expression values vs psi values for a 
#specific rMATs event
plot.result<-function(rbp,event,tissue)
{
  #get the ensembl ID for the target 
  events.tab<-read.table('/Users/karleg/Dimorph/fromGTF.SE.txt',header=T)

  event.label<-as.character(events.tab$GeneID[event])

  #read in the expression levels of RBPs
  rbp.exp<-read.table('/Users/karleg/Dimorph/rbp_expression0.txt',sep='\t',header=T)
  
  #keep the expression level of the relevant RBP only
  rbp.exp<-rbp.exp[grepl(rbp,rbp.exp$Name),]
  
  rbp.exp<-t(rbp.exp)
  
  colnames(rbp.exp)<-gsub("\\.[0-9]*","",rbp.exp[1,])
  
  rown<-rownames(rbp.exp)
  
  rbp.exp<-matrix(rbp.exp[-c(1,2),],ncol=1)
  
  rownames(rbp.exp)<-rown[-c(1,2)]
  
  rownames(rbp.exp)=gsub('\\.','-',rownames(rbp.exp))
  
  #read trhe GTEx meta data 
  meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)
  
  meta.data<-meta.data[meta.data$body_site_s %in% tissue,]
  
  rownames(meta.data)<-meta.data$Sample_Name_s
  
  #merge RBP expression and meta data tables
  rbp.exp<-merge(x=rbp.exp,y=meta.data,by='row.names')
  
  #Get the labels of the psi values which correspond to the GTEx individuals 
  psi.tab<-read.csv('/Users/karleg/Dimorph/rmats_final_matrix.se.jc.inc.txt',row.names = 1,header=TRUE,skip=event-1,nrows=1)
  
  psi.labels<-colnames(read.csv('/Users/karleg/Dimorph/rmats_final_matrix.se.jc.inc.txt',row.names = 1,header=TRUE,skip=0,nrows=1))
  
  psi.labels<-gsub("\\.","-",psi.labels)
  
  rbp.exp<-rbp.exp[rbp.exp$Row.names %in% psi.labels,]
  
  psi.labels.subset<-psi.labels[psi.labels %in% rbp.exp$Row.names]
  
  if (!(sum(rbp.exp$Row.names %in% psi.labels.subset)==nrow(rbp.exp) && sum( psi.labels.subset %in% rbp.exp$Row.names )==nrow(rbp.exp)))
    
    quit(status = 1)
  
  psi.vals<-as.numeric(psi.tab)
  
  psi.vals<-psi.vals[psi.labels %in% rbp.exp$Row.names]
  
  #order the psi values to match the RBP expression values
  psi.vals<-psi.vals[order(match(psi.labels.subset,rbp.exp$Row.names))]
  
  
  #Create two plots, one for male subjects and one for female subjects
  plot(as.numeric(as.character(rbp.exp[rbp.exp$sex_s=='male','V1'])),psi.vals[rbp.exp$sex_s=='male'],xlab='expression',ylab='psi',main=paste('Male',rbp,event.label))
  
  fit <- lm(psi.vals[rbp.exp$sex_s=='male']~as.numeric(as.character(rbp.exp[rbp.exp$sex_s=='male','V1'])))

  rmse <- round(sqrt(mean(resid(fit)^2)), 2)
  
  coefs <- coef(fit)
  
  b0 <- round(coefs[1], 4)
  
  b1 <- round(coefs[2],4)
  
  r2 <- round(summary(fit)$r.squared, 2)
  
  eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ r^2 == .(r2) * "," ~~ RMSE == .(rmse) ~~ Fvalp==.(1-pf(summary(fit)$fstatistic[1],df1=summary(fit)$fstatistic[2],df2=summary(fit)$fstatistic[3])))
  
  abline(fit)
  
  mtext(side=1,eqn) 
  
  plot(as.numeric(as.character(rbp.exp[rbp.exp$sex_s=='female','V1'])),psi.vals[rbp.exp$sex_s=='female'],xlab='expression',ylab='psi',main=paste('Female',rbp,event.label))
 
  fit <- lm(psi.vals[rbp.exp$sex_s=='female']~as.numeric(as.character(rbp.exp[rbp.exp$sex_s=='female','V1'])))
  
  rmse <- round(sqrt(mean(resid(fit)^2)), 2)
  
  coefs <- coef(fit)
  
  b0 <- round(coefs[1], 4)
  
  b1 <- round(coefs[2],4)
  
  r2 <- round(summary(fit)$r.squared, 2)
  
  eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ r^2 == .(r2) * "," ~~ RMSE == .(rmse)~~ Fvalp==.(1-pf(summary(fit)$fstatistic[1],df1=summary(fit)$fstatistic[2],df2=summary(fit)$fstatistic[3])))
  
  abline(fit)
  
  mtext(side=1,eqn) 
   
}

#create histograms of p-values

hist(unlist(level.tab))

hist(unlist(sex.tab))

#View the table of p-values
View(results[order(results[,1]),])

#plot RBP expression vs. arget tpsi for various RBP-target pairs
rbp='ENSG00000102081'

event=14263

tissue<-'Thyroid'

plot.result(rbp,event,tissue)

rbp='ENSG00000102081'

event=6267

plot.result(rbp,event,tissue)


rbp='ENSG00000125970'

event=1431

plot.result(rbp,event,tissue)

rbp='ENSG00000132819'

event=16639

plot.result(rbp,event,tissue)

rbp='ENSG00000102081'

event=27030

plot.result(rbp,event,tissue)

rbp='ENSG00000126945'

event=31657

plot.result(rbp,event,tissue)

#adipose

rbp='ENSG00000102317'

event=26648

plot.result(rbp,event,tissue)

rbp='ENSG00000132819'

event=26247

plot.result(rbp,event,tissue)


rbp='ENSG00000132819'

event=5735

plot.result(rbp,event,tissue)


rbp='ENSG00000102081'

event=13182

plot.result(rbp,event,tissue)




