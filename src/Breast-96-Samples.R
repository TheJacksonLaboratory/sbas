library(limma)
library(multtest)
library(Biobase)
library(edgeR)

#TODO:  Add arguments section to expose args to the user from command line
# Very nice R script example here: https://github.com/IARCbioinfo/needlestack/blob/8883a4ad65bb231f96f3e6e26095ec1b48be0842/bin/needlestack.r

# ############################## ARGUMENTS SECTION #############################
# ## Collect arguments
# args <- commandArgs(TRUE)
# 
# ## Parse arguments (we expect the form --arg=value)
# parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
# argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
# names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
# args <- argsL
# rm(argsL)
# 
# ## Give some value to options if not provided 
# if(is.null(args$opt_arg1)) {args$opt_arg1="default_option1"}
# if(is.null(args$opt_arg2)) {args$opt_arg2="default_option1"} else {args$opt_arg2=as.numeric(args$opt_arg2)}
# 
# ## Default setting when no all arguments passed or help needed
# if("--help" %in% args | is.null(args$arg1) | is.null(args$arg2)) {
#   cat("
#       The R Script arguments_section.R
#       
#       Mandatory arguments:
#       --arg1=type           - description
#       --arg2=type           - description
#       --help                - print this text
#       
#       Optionnal arguments:
#       --opt_arg1=String          - example:an absolute path, default:default_option1
#       --opt_arg2=Value           - example:a threshold, default:10
# 
#       WARNING : here put all the things the user has to know
#
#       Example:
#       ./arguments_section.R --arg1=~/Documents/ --arg2=10 --opt_arg2=8 \n\n")
#
#   q(save="no")
# }
#
# cat("first mandatory argument : ", args$arg1,"\n",sep="")
# cat("second mandatory argument : ", args$arg2,"\n",sep="")
# cat("first optional argument : ", args$opt_arg1,"\n",sep="")
# cat("second optional argument : ", args$opt_arg2,"\n",sep="")

# Thresholds
min_inclusion_threshold <- 0.25

# Body sites
body_site_of_choice <- 'Breast - Mammary Tissue'

# Samples to omit
toRemoveSampleNamePattern <- "11ILO"

#DS analysis
# "/Users/karleg/Downloads/rmats_final.se.jc.ijc.txt"
# "/Users/karleg/Downloads/rmats_final.se.jc.sjc.txt"
# "/Users/karleg/Downloads/SraRunTable.noCram.noExome.noWGS.totalRNA.txt"

# Input matrices, assuming we work in pairs of inclusion/skipping exons
inclusion_matrix_path <- "testdata/rmats_final.se.jc.ijc.txt.gz"
skipping_matrix_path  <- "testdata/rmats_final.se.jc.sjc.txt.gz"
meta_data             <- "testdata/SraRunTable.noCram.noExome.noWGS.totalRNA.txt"

# SRA ids are the colnames 
inc.iso.counts.mem  <- readr::read_csv(inclusion_matrix_path)
skip.iso.counts.mem <- readr::read_csv(skipping_matrix_path)

# SRA ids are the entries in the column named 'Run', 
meta.data           <- readr::read_csv(meta_data)
# Keep only SRA samples in metadata that are also present in the inclusion rMATS matrix
meta.data     <- meta.data[meta.data$Run %in% colnames(inc.iso.counts.mem),]

dimensionsMatch <- function() {
  n_sample_inc   <- ncol(inc.iso.counts.mem [, !grepl("ID", colnames(inc.iso.counts.mem))])
  n_sample_skip  <- ncol(skip.iso.counts.mem[, !grepl("ID", colnames(skip.iso.counts.mem))])
  n_sample_meta  <- nrow(meta.data)
  
  if (!(all(sapply(list(n_sample_inc,n_sample_skip,n_sample_meta), function(x) x == n_sample_meta)))) {
    return(FALSE)
    stop("Samples across rMATS matrices and metatada is not the same")
  }
  if ((all(sapply(list(n_sample_inc,n_sample_skip,n_sample_meta), function(x) x == n_sample_meta)))) {
    message("All checked! Samples across rMATS matrices and metatada match")
    return(TRUE)
  }
}


# CHECKPOINT: 
# Check if rMATS and metadata dimensions match; if not stop and provide error message

# Check if samples names in metadata and the 2 rMATS matrices match
if ( !(isTRUE(dimensionsMatch())) ) {
  stop("Samples across rMATS matrics and metatada is not the same")
}
if ( (isTRUE(dimensionsMatch())) ) {
  message("\nApplying samples filtering to rMATS and metadata tables ..\n")
}

# END OF CHECKPOINT


#_______ inc ____________

# Name rMATS matrix rows and remove ID column by name
rownames(inc.iso.counts.mem) <- inc.iso.counts.mem$ID
inc.iso.counts.mem$ID <- NULL

# Order rMATS matrix based on metadata order
inc.iso.counts.mem <- inc.iso.counts.mem[, order(match(colnames(inc.iso.counts.mem), meta.data$Run))]

# Filter samples in rMATS matrix:
# 1. remove exclusion list samples 'toRemoveSampleNamePattern'
# 2. keep samples that correspond to tissue under examination 'body_site_of_choice'
toKeepSraIds <- colnames(skip.iso.counts.mem)[!grepl(toRemoveSampleNamePattern, meta.data$Sample.Name) & (meta.data$body_site %in% body_site_of_choice)]
inc.iso.counts.mem<-inc.iso.counts.mem[ , toKeepSraIds]


#_______ skip ____________
rownames(skip.iso.counts.mem) <- skip.iso.counts.mem$ID
skip.iso.counts.mem$ID <- NULL

# Order rMATS matrix based on metadata order
skip.iso.counts.mem <-skip.iso.counts.mem[,order(match(colnames(skip.iso.counts.mem), meta.data$Run))]

# Filter samples in rMATS matrix:
# 1. remove exclusion list samples 'toRemoveSampleNamePattern'
# 2. keep samples that correspond to tissue under examination 'body_site_of_choice'
skip.iso.counts.mem <-skip.iso.counts.mem[,toKeepSraIds]

#_______ metadata ____________

# Remove sample ids that match the exclusion list 'toRemoveSampleNamePattern'
meta.data<-meta.data[!grepl(toRemoveSampleNamePattern,meta.data$Sample.Name),]

# Keep sample ids that match the exclusion list 'toRemoveSampleNamePattern'
meta.data<-meta.data[meta.data$body_site %in% body_site_of_choice,]

# Keep samples with sum of counts across features
keep.cols <- colSums(inc.iso.counts.mem > 0) >= min_inclusion_threshold * nrow(inc.iso.counts.mem)
keep.cols <- keep.cols & colSums(skip.iso.counts.mem > 0) >= min_inclusion_threshold*nrow(skip.iso.counts.mem)

inc.iso.counts.mem=inc.iso.counts.mem[,keep.cols]
skip.iso.counts.mem=skip.iso.counts.mem[,keep.cols]
meta.data=meta.data[keep.cols,]

if (sum(colnames(skip.iso.counts.mem)==meta.data$Run)!=nrow(meta.data))
   stop("Metadata dimensions don't match input data dimensions. Check: sum(colnames(skip.iso.counts.mem)==meta.data$Run)!=nrow(meta.data)")

if (sum(colnames(inc.iso.counts.mem)==meta.data$Run)!=nrow(meta.data))
  stop("Metadata dimensions don't match input data dimensions. Check: sum(colnames(inc.iso.counts.mem)==meta.data$Run)!=nrow(meta.data)")

exprs.mat <- matrix(data=cbind(as.matrix(skip.iso.counts.mem),
                               as.matrix(inc.iso.counts.mem)),
                    ncol=2*ncol(skip.iso.counts.mem),
                    nrow=nrow(skip.iso.counts.mem),
                    dimnames=list(1:nrow(skip.iso.counts.mem),
                                  c(colnames(skip.iso.counts.mem),
                                    paste(colnames(skip.iso.counts.mem),'2',sep=''))))

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

write.table(res[res$adj.P.Val<=0.05 & abs(res$logFC)>=log2(1.5),],paste(paste(body_site_of_choice,collapse='.'),'se.txt',sep=''),row.names = T,col.names = T,quote = F)

write.table(rownames(res),paste0(paste(body_site_of_choice,collapse='.'),'all_genes.txt'))


#Run Ontologizer

meta.data<-read.table('/Users/karleg/Downloads/fromGTF.SE.txt',sep='\t',header=TRUE)
de.tab<-read.table('Breast - Mammary Tissuese.txt')
de.tab.with.meta<-merge(de.tab,meta.data,by.x='row.names',by.y='ID')
all.genes<-read.table('Breast - Mammary Tissueall_genes.txt')
write.table(meta.data$geneSymbol[meta.data$ID %in% all.genes[,1]],'universe.txt',quote = F,row.names = F,col.names = F)
write.table(de.tab.with.meta$geneSymbol,'gene_set.txt',quote = F,row.names = F,col.names = F)
system('java -jar /Users/karleg/Ontologizer/Ontologizer.jar -g /Users/karleg/Ontologizer/go.obo -a /Users/karleg/Ontologizer/goa_human.gaf -s gene_set.txt -p universe.txt -c Term-For-Term -m Benjamini-Hochberg -n')


# unique(meta.data$body_site)
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
# > dim(meta.data)
# [1] 185  79