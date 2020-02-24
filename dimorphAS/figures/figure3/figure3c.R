## This script will produce Figure 3c
## The first part of the script was used to generate the analysis
## To create the figure, skip down and read in the datafile.


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

tissue <- tissue.list[[2]]

load(paste('/Users/karleg/Dimorph/McmcMostVaryingMoreSigs_',tissue,'.Rdata',sep=''))

mcmcCoda<-mcmcCoda[,which(grepl('beta3\\[101\\]',varnames(mcmcCoda))),drop=FALSE]

diagMCMC( mcmcCoda , parName=c("beta3[101]") )  

save.image('/Users/karleg/Dimorph/RDATA/figureS3.RData')

## Create Figure 3c -- genes with >10 AS events in multiple tissues
# Removal all variables from workspace
rm(list=ls())
#Before running the following, use the Session menu to set working directory to source file location
#This file will save a PDF called hdiExample.pdf
load('figure3c.RData')
library(coda)

npgBlue<- rgb(60/256,84/256,136/256,1)
npgRed <- rgb(220/256,0,0,0.8)
npgGreen <- rgb(0,160/256,135/256,1)
npgBrown <- rgb(126/256,97/256,72/256,1)


myDbdaDensPlot = function( codaObject , parName=varnames(codaObject)[1] , plColors=NULL ) {
  if ( all( parName != varnames(codaObject) ) ) { 
    stop("parName must be a column name of coda object")
  }
  nChain = length(codaObject) # or nchain(codaObject)
  if ( is.null(plColors) ) plColors=1:nChain
  xMat = NULL
  yMat = NULL
  hdiLims = NULL
  for ( cIdx in 1:nChain ) {
    densInfo = density(codaObject[,c(parName)][[cIdx]]) 
    xMat = cbind(xMat,densInfo$x)
    yMat = cbind(yMat,densInfo$y)
    hdiLims = cbind(hdiLims,HDIofMCMC(codaObject[,c(parName)][[cIdx]]))
  }
  matplot( xMat , yMat , type="l" , col=plColors , 
           main="" , xlab="Param. Value" , ylab="Density" )
  abline(h=0)
  points( hdiLims[1,] , rep(0,nChain) , col=plColors , pch="|" )
  points( hdiLims[2,] , rep(0,nChain) , col=plColors , pch="|" )
  text( mean(hdiLims) , 0 , "95% HDI" , adj=c(0.5,-0.2), cex=2 )
  EffChnLngth = effectiveSize(codaObject[,c(parName)])
  MCSE = sd(as.matrix(codaObject[,c(parName)]))/sqrt(EffChnLngth) 
  text( max(xMat) , max(yMat) , adj=c(1.0,1.0) , cex=1.75 ,
        paste("MCSE =\n",signif(MCSE,3)) )
}

mydiag <- function( codaObject , parName=varnames(codaObject)[1] ,
                    saveName=NULL , saveType="pdf" ) {
  DBDAplColors =  c("skyblue","black","royalblue","steelblue")#c(npgBlue,"black",npgRed,npgGreen)
  openGraph(height=5,width=7)
  par( mar=0.5+c(3,4,1,0) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
       cex.lab=2 , cex.axis=2)
  require(coda)
  # The following is the HDI plot
  myDbdaDensPlot(codaObject,parName,plColors=DBDAplColors)
  #Uncomment the following to add the title (beta3[101]), which is not informative so
  # we will leave it out.
  #mtext( text=parName , outer=TRUE , adj=c(0.5,0.5) , cex=2.0 )
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName), type=saveType)
  }
}


mydiag( mcmcCoda , parName=c("beta3[101]") , saveName="fig3c")  

