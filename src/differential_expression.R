#!/usr/bin/env Rscript

############################## ARGUMENTS SECTION #############################
## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no all arguments passed or help needed
if("--help" %in% args | "help" %in% args | (length(args) == 0) | (length(args) == 1) ) {
  cat("
      The R Script differential_expression.R

      Mandatory arguments:
          --counts=path             - A comma seperated file with the counts from rMATS. 
                                      Rownames must be the features (e.g. splicing events)
                                      Column names must be the 'Run' ids from the 'metadata' file (e.g. 'SRR1068788'). 
                                      The ID column must be present in the file.
                                      Compressed .gz format accepted (e.g. rmats_final.se.jc.ijc.txt.gz)
                                      
          --metadata=path           - A comma seperated file with metadata about the samples in the 'counts' file.
                                      Rownames must be the samples.
                                      Column names must be metadata about the samples (e.g. 'Sample.Name','Run', 'body_site', etc). 
                                      The 'Sample.Name' column must be present in the file and the values of the columns
                                      must include the colnames() in the 'counts' file.
                                      Compressed .gz format accepted (e.g. SraRunTable.noCram.noExome.noWGS.totalRNA.txt.gz)

          --tissue=value            - A string in single quotes. Must be a value from the column 'body_site' in the metadata file. 
                                      (e.g. 'Breast - Mammary Tissue')

         --help                     - you are reading it

      Optionnal arguments:
          --threshold=Value         - a threshold, default:0.25
          
          --exclusion_list=path     - A single column file without a header with the 'Run' ids (e.g. SRR1361860) from the 'metadata' file 
                                      to exclude from the analysis. An example with the 'Run' ids that correspond to the sample with 
                                     'biospecimen repository sample id' can be found in 'sbas/testdata/exlude.txt'
                                      Compressed .gz format accepted (e.g. SraRunTable.noCram.noExome.noWGS.totalRNA.txt.gz)

      Usage:
      
          The typical command for running the script is as follows:
    
          ./differential_expression.R --counts='testdata/counts.csv' --metadata='testdata/metadata.csv' --threshold=0.20 
      
      
      WARNING : here put all the things the user has to know
      \n")
  
  
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")

argsL        <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args         <- argsL
rm(argsL)

## Give some value to optional arguments if not provided
if(is.null(args$threshold)) {args$threshold =0.25} else {args$threshold=as.numeric(args$threshold)}

cat("\n")
cat("ARGUMENTS SUMMARY")
cat("\n")
cat("counts          : ", args$counts,  "\n",sep="")
cat("metadata        : ", args$metadata, "\n",sep="")
cat("tissue          : ", args$tissue,   "\n",sep="")
cat("threshold       : ", args$threshold,"\n",sep="")
if(!(is.null(args$exclusion_list))) { cat("exclusion_list  : ", args$exclusion_list,"\n",sep="") }
if(!(is.null(args$threshold))) {args$threshold =0.25} else {args$threshold=as.numeric(args$threshold)}

############################## LIBRARIES SECTION #############################

suppressWarnings(suppressMessages(library(limma)))
suppressWarnings(suppressMessages(library(multtest)))
suppressWarnings(suppressMessages(library(Biobase)))
suppressWarnings(suppressMessages(library(edgeR)))

############################### SCRIPT SECTION ###############################


if(!(is.null(args$exclusion_list))) {
  exclusion_list <- as.vector(t(readr::read_csv(file = args$exclusion_list, )))
}
if(is.null(args$exclusion_list)) {
  exclusion_list <- ''
}

counts   <- readr::read_csv(file = args$counts)
metadata <- readr::read_csv(file = args$metadata)
tissue   <- as.character(args$tissue)


# TODO remove this
exclusion_list <- as.vector(t(readr::read_csv(file = "testdata/exlude.txt", col_names = FALSE)))
counts   <- readr::read_csv(file = "testdata/rmats_final.se.jc.ijc.txt.gz")
metadata <- readr::read_csv(file = "testdata/SraRunTable.noCram.noExome.noWGS.totalRNA.txt")
tissue   <- 'Breast - Mammary Tissue'
dim(counts)

# Select SRA ids that correspond to specified --tissue
toKeepSraIdsTissue <- metadata$Run[ metadata$body_site == tissue]
length(toKeepSraIdsTissue)

# Filter out further exclusion list SRA ids if specified
toKeepSraIds <- toKeepSraIdsTissue[!(toKeepSraIdsTissue %in% exclusion_list)]
length(toKeepSraIds)

# Select all toKeepSraIds included in the counts matrix (not all samples in metadata are available)
toKeep <- colnames(counts)[colnames(counts)  %in% toKeepSraIds]
toKeep <- c('ID', toKeep)
length(toKeep)

# TODO elliminate printing above
# Subset counts, metadata based on all exclusion, inclusion criteria 
counts   <-  as.data.frame( counts[ , toKeep ])
metadata <- as.data.frame(metadata[metadata$Run  %in% toKeep, ])


# Name rMATS matrix rows and remove ID column by name
rownames(counts) <- counts$ID
counts$ID <- NULL

# Order rMATS matrix based on metadata order
counts <- counts[, order(match(colnames(counts), metadata$Run))]


exprs.mat <- matrix(data=cbind(as.matrix(skip.iso.counts.mem),
                               as.matrix(counts)),
                    ncol=2*ncol(skip.iso.counts.mem),
                    nrow=nrow(skip.iso.counts.mem),
                    dimnames=list(1:nrow(skip.iso.counts.mem),
                                  c(colnames(skip.iso.counts.mem),
                                    paste(colnames(skip.iso.counts.mem),'2',sep=''))))

keep.events<-rep(T,nrow(exprs.mat))

groups=c(paste0(metadata$sex,'-skip'),paste0(metadata$sex,'-inc'))

groups<-groups[colSums(exprs.mat)>0]

exprs.mat<-exprs.mat[,colSums(exprs.mat)>0]

zero.rows=rep(F,nrow(exprs.mat))

for (group in c('male-skip','female-skip','male-inc','female-inc'))

  zero.rows=zero.rows | (rowSums(exprs.mat[,groups %in% group] >0) < (1/4)*min(table(groups)))

exprs.mat[zero.rows,]=0

for (group in c('male-skip','female-skip','male-inc','female-inc'))

  keep.events<-keep.events & (rowSums((exprs.mat[,groups %in% group]) >= 1) >= (1/4)*min(table(groups)))

rm(exprs.mat)

inc.counts.male<-counts[,metadata$sex=='male']

inc.counts.female<-counts[,metadata$sex=='female']

skip.counts.male<-skip.iso.counts.mem[,metadata$sex=='male']

skip.counts.female<-skip.iso.counts.mem[,metadata$sex=='female']

inc.counts.mat<-cbind(inc.counts.male,inc.counts.female)

rownames(inc.counts.mat)<-rownames(counts)

skip.counts.mat<-cbind(skip.counts.male,skip.counts.female)

rownames(skip.counts.mat)<-rownames(skip.iso.counts.mem)

counts.mat<-cbind(inc.counts.male,inc.counts.female,skip.counts.male,skip.counts.female)

rownames(counts.mat)<-rownames(counts)

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

write.table(res[res$adj.P.Val<=0.05 & abs(res$logFC)>=log2(1.5),],paste(paste(body_site_of_choice,collapse='.'),'se.txt',sep=''),row.names = T,col.names = T,quote = F)

write.table(rownames(res),paste0(paste(body_site_of_choice,collapse='.'),'all_genes.txt'))


#Run Ontologizer

metadata<-read.table('/Users/karleg/Downloads/fromGTF.SE.txt',sep='\t',header=TRUE)
de.tab<-read.table('Breast - Mammary Tissuese.txt')
de.tab.with.meta<-merge(de.tab,metadata,by.x='row.names',by.y='ID')
all.genes<-read.table('Breast - Mammary Tissueall_genes.txt')
write.table(metadata$geneSymbol[metadata$ID %in% all.genes[,1]],'universe.txt',quote = F,row.names = F,col.names = F)
write.table(de.tab.with.meta$geneSymbol,'gene_set.txt',quote = F,row.names = F,col.names = F)
system('java -jar /Users/karleg/Ontologizer/Ontologizer.jar -g /Users/karleg/Ontologizer/go.obo -a /Users/karleg/Ontologizer/goa_human.gaf -s gene_set.txt -p universe.txt -c Term-For-Term -m Benjamini-Hochberg -n')


# unique(metadata$body_site)
# [1] Breast - Mammary Tissue
# 55 Levels:
# > skip.iso.counts.mem[1:3,1:3]
# SRR821498 SRR808428 SRR808942
# 1         0         1         0
# 2         0         0         1
# 3         0         0         0
# > inc.iso.counts.mem[1:3,1:3]
# SRR821498 SRR808428 SRR808942
# 1         0         0         0
# 2       203       202       396
# 3         0         2         0
# > dim(exprs.mat)
# [1]  39 370
# > dim(inc.iso.counts.mem)
# [1]  39 185
# > dim(skip.iso.counts.mem)
# [1]  39 185
# > dim(metadata)
# [1] 185  79