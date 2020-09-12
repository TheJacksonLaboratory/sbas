# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# # differenticalSplicingJunctionExpressionAnalysis
#
# rMATS 3.2.5 was run on controlled access RNASeq files retrieved experiments stored in the Sequence Read Archive with controlled access managed by dbGaP.   This experiment was run with the fastq files from GTEx v8.
#
# The output (read in by section 1.2.1) are matrices which are the result of executing the rmats-nf nextflow workflow on all the samples from GTEx V8 https://github.com/lifebit-ai/rmats-nf.   The workflow begins with the accessions file and continues until a matrix.  Run without statistics, for the purposes of rMATS creating an annotated junction file for each of the five (5) splicing types.  The matrix is possible with this version of rMATS, as the junction ID is unique per annotation GTF.  In this running, we used gencode.v30.annotation.gtf (complete annotation).   The result is 5 matrices per splicing type.
#
# rMATS RNASeq-MATS.py produces 10 different output types which get assembled into as type junction ID by sample ID matrices
#
# ## Alternative Splice Site Types are: (se, a3ss, a5ss, mxe, ri)
#
#   * Skipped Exon events (se),
#   * Alternative 3' splice site (a3ss),
#   * Alternative 5' splice site (a5ss),
#   * Mutually exclusive exon (mxe),
#   * and retention intron (ri)
#
# ## There are two different kinds of junction counts
#
# For our analysis here, we used just the jc count matrices.
#   * jc = junction counts - reads that cross the junction
#   * jcec = junction counts plus reads on the target (such as included exon)
#
# ## And the count type -- there are 5 types
#
#   * inclusion levels (percent spliced in)
#   * included junction counts (ijc)
#   * skipped junction counts (sjc)
#   * inclusion length (inclen)
#   * skipped length (skiplen)
#
# ## 1. Loading dependencies

suppressMessages({
    options(warn = -1) 
    library(readr)
    library(edgeR)
    library(limma)
    library(multtest)
    library(Biobase)
    library(tibble)
    library(R.utils)
    library(snakecase)
})

# ## 1.1 Nextflow execution parameter execution
#
# Using the papermill library, we can parallelize execution of this notebook.  To do this, the loop at the bottom of this notebook should be commented out -- and papermill will run across all tissues fed into it.
# parameters for nextflow execution of notebook
# tissue_index = 17

# ### 1.2 get rMATS GTF annotations
#
# For each splicing type, the junctions are defined, so we have 5 specific annotated splicing specific junction ID annotation files:
# fromGTF.A3SS.txt <- annotations for the alternative 3' splice site junctions
# fromGTF.A5SS.txt <- annotations for the alternative 5' splice site junctions
# fromGTF.MXE.txt <- annotations for the mutually exclusive exon junctions
# fromGTF.RI.txt <- annotations for the retained introns junctions
# fromGTF.SE.txt <- annotations for the skipped exon junctions

getReleasedGTFAnnotations <- function ( destDir ) {
    if (file.exists("../data/fromGTF.A3SS.txt")) {
        message("Found unpacked fromGTF.* files")
    } else {
        system("mkdir -p ../data && tar xvfz ../data/fromGTF.tar.gz -C ../data", intern = TRUE)
        message("Done Decompressing fromGTF.tar.gz into ../data")
        system("gunzip ../data/fromGTF.*txt.gz", intern = TRUE)
        message("Done Gunzipping files into ../data")
    }
}

# ## 1.2.3 Read in SraRunData metadata 
#
# - `Sequence Read Archive (SRA)` Accession Data, `SRR` numbers, this is used to map the SRR accession numbers to the sample information (SAMPID) which will be used to obtain the phenotype information.

getSraRunData <- function ( destDir ) {
  message("Loading metadata from SraRunTable.txt.gz ../data/ ..\n")
  srr_metadata <- data.table::fread("../data/SraRunTable.txt.gz")
  srr_metadata$SAMPID <- gsub('-','\\.',srr_metadata$biospecimen_repository_sample_id)
  message("done!\n")
  return(srr_metadata)
}
# ### 1.2.4 Renew GTEx expression object
#
# - `Genome Tissue Expression (GTEx)` Clinical Annotation - this is the expressionSet object that has the phenotype information
# In this analysis there will be 3 expressionSet objects.  This one contains the gene Expression Count Data and phenotypes  -- this is done in two parts, pull it from yarn, and then inspect and correct.

inspectAndCorrectExpressionSetObject <- function ( es ) {
    # This function checks whether the samples represented by
    # colnames(exprs(es) are the same as rownames(pData(es))
    # if not, we take the intersect
   
   sample_names=as.vector(as.character(colnames(exprs(es))))
   pheno_sample_names=as.vector(as.character(rownames(pData(es))))
   intersecting_samplenames <- intersect(sample_names, pheno_sample_names)

   message("sample_names n=", length(sample_names), " pheno_sample n=", length(pheno_sample_names),
           "intersecting_samplenames n=", length(intersecting_samplenames))

   logical_match_names=superset %in% subset
   logical_match_names <- rownames(pData(es)) %in% intersecting_samplenames
   message("\nLogical diff pheno_sample_names, expression_sample_names\n",
        paste(table(logical_match_names)), collapse = " ")
#######
   pData(es) <- pData(es)[logical_match_names==TRUE,]

   message("\nAFTER: dim(exprs(es))[2]\n",
        paste(dim(exprs(es))[2]), collapse=" ")
   message("\nAFTER:: dim(pData(es))[1]\n",
        paste(dim(pData(es))[1]), collapse=" ")
   message("\nAFTER:: dimension of es\n",
        paste(dim(es)[2]), collapse=" ")
   message("Saving corrected gtex.rds..\n")
   saveRDS(es, file="../data/gtex.corrected.rds")   
   message("Generating sha256sum for gtex.corrected.rds ..\n")    
   message(system("sha256sum ../data/gtex.corrected.rds", intern = TRUE))
   message("Done!\n")

   return(es)
}

# +
renewGTExExpressionSet <- function ( destDir ) {
  if (!("gtex.rds" %in% list.files(destDir))) {
      message("Downloading and loading obj with GTEx v8 with 'yarn::downloadGTExV8()'\n")
      es <- yarn::downloadGTExV8(type='genes',file='../data/gtex.rds')
      message("Done!\n")
      message("Generating sha256sum for gtex.rds ..\n")    
      message(system("sha256sum ../data/gtex.rds", intern = TRUE))
      message("Done!\n")
      message("running inspectAndCorrectExpressionSet to generate gtex.corrected.rds\n")
      es <- inspectAndCorrectExpressionSetObject (es)
      message("Done!\n")
      
  } else {
    # Load with readRDS() if gtex.rds available in data/
      message("Loading obj GTEx v8 corrected rds object with readRDS from ../data/gtex.corrected.rds\n")   
      es <- readRDS(file = "../data/gtex.corrected.rds")
      message("Done!\n")
      message("Generating sha256sum for gtex.corrected.rds ..\n")    
      message(system("sha256sum ../data/gtex.corrected.rds", intern = TRUE))
      message("Done!\n")
  }
  return (es)

}
# -

# ### 1.2.5 get gtex.corrected.rds

# +
getGTExCorrectedRDS <- function ( destDir ) {

  message("Loading gtex.corrected.rds ..\n")
  es <- readRDS(file = "../data/gtex.corrected.rds")
  message("done!\n")
  return(es)

}
# -

# ### 1.2.6 get reduced Tissue Data
#
# Stored in the assets subdirectory, reduced by inspection and selection focusing on those tissues with sufficient samples.

getTissueReduction <- function ( filename ) {

    tissue_reduction <- read.table(filename, header=TRUE, sep="\t",
                               skipNul=FALSE, stringsAsFactors = FALSE)
    colnames(tissue_reduction)  <- c("SMTSD","female","male","include","display_name")

    return(tissue_reduction)
}
# ### 1.2.7 make phenotype data matrix portion of ExpressionSet (makePData)
#
# Matching SraRunTable, Sequence Read Archive Sequence Run accession number, SRR to the sample, SAMPID, that will be used to match the SRR to the samples phenotype information to make the phenotype data portion of the expressionSet object.
#
# counts data supplies the SRR through the column
# Metadata from the SraRunTable.txt file matches the SRR to the SAMPID
# Phenotype data retrieved from the GTEx gene ExpressionSet data object for that SAMPID
#
# SRR -> SAMPID -> Phenotype Data
#
# This routine takes a while to run, could be improved by changing it to a lapply function.
# Luckily, this needs only to be run once for all 10 matrices since the ijc, sjc counts were done on the same samples and defined alternative splicing events.
makePData <- function (counts_srr, gtexPhenoDataObj, srr_metadata) {
    started = FALSE
    srr_missing = FALSE
    pd_missing = FALSE
    for (i in 1:length(counts_srr)) {
        #
        # SraRunTable connects SRR to SAMPID
        #
        srr_metadata_match <- as.character(srr_metadata$Run) %in% as.character(counts_srr[i])
        if (sum(srr_metadata_match) < 1){
            srr_metadata_missing <- as.character(srr_metadata$Run) 
            message("no match srr_metadata$Run! ERROR!")
            if (!srr_missing) {
                srr_metadata_missing_list <- srr_metadata_missing
                srr_missing = TRUE
            } else {
                srr_metadata_missing_list <- rbind(srr_metadata_missing_list, srr_metadata_missing)
            }
        } else {
            srr_metadata_sampid <- srr_metadata[srr_metadata_match==TRUE,]$SAMPID
#            message("\nsrr_metadata_sampid\n",
#               paste(srr_metadata_sampid), collapse=" ")
            #
            # SAMPID is used to match SRR to Phenotype Data for GTEx
            #
            pd_sampid_match <- as.character(pData(gtexPhenoDataObj)$SAMPID) %in% as.character(srr_metadata_sampid)
#           message("\npd_sampid_match table\n",
#               paste(table(pd_sampid_match)), collapse=" ")
#           message("\nsum(pd_sampid_match)\n",
#               paste(sum(pd_sampid_match)), collapse = " ")
            if (sum(pd_sampid_match) < 1 ) {
#               message("no pd_sampid match to srr_metadata_sampid!\n")
                pd_sampid_missing <- as.character(srr_metadata_sampid)
                if (!pd_missing) {
                    pd_sampid_missing_list <- pd_sampid_missing
                    pd_missing = TRUE
                } else {
                    pd_sampid_missing_list <- rbind(pd_sampid_missing_list, pd_sampid_missing)
                }
            }  else {
                pd     <- pData(gtexPhenoDataObj)[pd_sampid_match,]
                srr_metadata_srr <- srr_metadata[srr_metadata_match==TRUE,]$Run
                pd$SRR <- srr_metadata_srr
                if (started == TRUE) {
                    pdfinal <- rbind(pdfinal, pd)
                } else {
                    pdfinal <- pd
                    started = TRUE
                }
            } # end if pd_sampid_match
        } # end if srr_metadata_match
    }# end for loop

    pdfinal_filename           <- paste0(paste0("../data/", counts_name), "pdata.csv")
    if (srr_missing) {
        metadata_missing_filename  <- paste0(paste0("../data/", counts_name), "srr_metadata_missing.csv")
        write.csv(srr_metadata_missing_list, metadata_missing_filename , quote=FALSE, row.names=FALSE)
    }
    if (pd_missing) {
        pd_sampid_missing_filename <- paste0(paste0("../data/", counts_name), "pd_sampid_missing.csv")
        write.csv(pd_sampid_missing_list   , pd_sampid_missing_filename, quote=FALSE, row.names=FALSE)
    }
    write.csv(pdfinal                  , pdfinal_filename          , quote=FALSE, row.names=FALSE)
    return(pdfinal)
}

# ### 1.2.8 get GTEx phenotype data for the SRR accessions
# Transitive closure permits the association of the sequence reads, SRR Accessions, through the SraRunTable.txt (obtained from selecting annotation from the dbGaP login) with the SAMPID used with the GTEx.  This SAMPID is the means by which we can get this phenotype data and associate it with the counts data.

getGTExPhenoDataForSRR <- function (destDir) {

    message("Loading srr_pdata\n")
    srr_pdata <- readr::read_csv("../data/srr_pdata.csv")
    message("done!\n")

    return(srr_pdata)
}
# ### 1.2.9 makeCountsMatrix 
#
# Given the counts filename, make a data matrix.

makeCountsMatrix <- function (filename_gz) {
    message("\nloading ", paste(filename_gz), collapse=" ")
    counts <- data.table::fread(filename_gz)
    message("done!")
    rownames(counts) <- counts$ID
    counts <- counts[,-1]
    counts <- data.matrix(counts)
    return(counts)
}

#
# ### 1.2.10 makeSplicingExpressionSetObject 
#
# Given the phenotype data object for each of the count matrices, create an expressionSet object for each to facilitate analysis.

makeSplicingExpressionSetObject <- function (srr_pdata, counts) {
    message("making splicing expressionSet object")
    #
    # match srr counts with srr_pdata - there were some srr without phenodata
    #
    pdata_match <- as.character(colnames(counts)) %in% as.character(srr_pdata$'SRR')
    
    counts <- counts[,pdata_match]

    #
    # reorder the srr_pdata to match the colnames of the counts
    #
    reorder_idx <- match(as.character(colnames(counts)), as.character(srr_pdata$'SRR'))
    srr_pdata   <- srr_pdata[reorder_idx,]

    #
    # make the srr_pdata an AnnotatedDataFrame
    #
    metadata <- data.frame(labelDescriptions=as.character(colnames(srr_pdata)))
    phenoData <- new("AnnotatedDataFrame", data = srr_pdata, varMetadata=metadata)

    #
    # make the counts an expressionSet
    # and provide the phenoData (the annotatedDataFrame construct from above)
    #
    es <- ExpressionSet(as.matrix(counts))
    phenoData(es) <- phenoData
    
    message("made new expressionSet object\n",
           paste(dim(es)), collapse = " ")
    message("dim pData(es)\n",
           paste(dim(pData(es))), collapse = " ")
    message("dim exprs(es)\n",
           paste(dim(exprs(es))), collapse = " ")
    message("done!\n")
    
    return(es)
}

# ## 1.3 Preprocessing 
#

# ### 1.3.1  Reduce Sample Set 
# Read in all requirements so that the stage is properly set -- tissues.tsv contains the subset of files desired for analysis.
# It is found in the `assets` subdirectory

reduceSampleSet <- function (tissue_reduction, es) {

   message("\nsize tissue_reduction\n",
        paste(dim(tissue_reduction), collapse=" "))
   message("\nsize es\n",
        paste(dim(es)), collapse=" ")
   message("\nsize pData(es)\n",
        paste(dim(pData(es)), collapse=" "))
   # only include those tissues we wish to continue with
   message("\n number of tissue types to keep\n",
        paste(table(tissue_reduction$include)), collapse = " ")
   tissue_reduction <- tissue_reduction[tissue_reduction$include==1,]

   # create a matching tissue name to go with the expressionSet phenotype esect
   pData(es)$SMTSD        <- factor(snakecase::to_snake_case(as.character(pData(es)$SMTSD)))
   tissue_reduction$SMTSD <- factor(snakecase::to_snake_case(as.character(tissue_reduction$SMTSD)))

   message("\nlength tissues in phenotype data\n",
        paste(length(levels(pData(es)$SMTSD)), collapse = " "))
   message("\nlength tissues in tissue_reduction data\n",
        paste(length(tissue_reduction$SMTSD), collapse = " "))

   keep <- pData(es)$SMTSD %in% tissue_reduction$SMTSD
   message("\nhow many to keep in phenotype data\n",
        paste(table(keep), collapse = " "))

   es        <- es       [          ,keep==TRUE]
   pData(es)$SMTSD <- factor(pData(es)$SMTSD)
   message("\nsize reduced es\n",
        paste(dim(es)), collapse=" ")
   message("\nsize pData(es)\n",
        paste(dim(pData(es)), collapse=" "))
   message("\nsize exprs(es)\n",
        paste(dim(exprs(es)), collapse = " "))

   # test to make sure we don't have nonsense
   keep = pData(es)$SMTSD== "breast_mammary_tissue"
   message("\nTEST: how many to keep in to have only breast_mammary_tissue\n",
        paste(table(keep), collapse = " "))
   tes        = es       [          ,keep==TRUE]
   pData(tes) = pData(es)[keep==TRUE,          ]
   message("\nTEST: size breast_mammary_tissue tes\n",
        paste(dim(tes), collapse=" "))
   message("\nTEST: size phenotype object pData(tes)\n",
        paste(dim(pData(tes)), collapse=" "))
   # end test
   return (es)
}


eliminateChrYfromGTF <- function ( fromGTF ) {

   fromGTF.keepAllButChrY <- (fromGTF$chr != "chrY")
   fromGTF           <- fromGTF[fromGTF.keepAllButChrY,]
   rownames(fromGTF) <- fromGTF$ID
   return(fromGTF)
}

eliminateChrYfromExpressionSet <- function ( fromGTF, es ) {

    es_row_ids <- rownames(es)
    gtf_ids <- fromGTF$ID
    keep <- es_row_ids %in% gtf_ids
   
    es   <- es[keep == TRUE, ]
    message("remade the expressionSet object\n",
           paste(dim(es)), collapse = " ")
    message("dim pData(es)\n",
           paste(dim(pData(es))), collapse = " ")
    message("dim exprs(es)\n",
           paste(dim(exprs(es))), collapse = " ")
    message("done!\n")

    return(es)
}

# ### 2.1 Differential analysis as_event (combined ijc and sjc)
#
# Differential Analysis (DE) was performed using voom (Law et.al., 2014) to transform junction counts (reads that were aligned to junctions when an exon is included - ijc, and reads that were aligned to junctions when the exon is excluded - sjc) with associated precision weights, followed by linear modeling and empirical Bayes procedure using limma.    In each tissue, the following linear regression model was used to detec secually dimorphic alternative splicing event expression: 
#
#            y = B0 + B1 sex + B2 as_event + B3 sex*as_event + epsilon (error)  
#            
#
# where y is the alternative splicing event expression; sex denotes the reported sex of the subject, as_event represents the specific alternative splicing event - either included exon junction counts or skipped exon junction counts and their interaction terms.   Donor is added to our model as a blocking variable used in both the calculation of duplicate correlation as well as in the linear fit.

# ### 2.2 Voom, limma's lmFit and eBayes
#
# Using sample as a blocking variable, we are able to model the effects of the donor on the results, which improves the power.  This topic is discussed in biostars https://www.biostars.org/p/54565/.  And Gordon Smyth answers the question here https://mailman.stat.ethz.ch/pipermail/bioconductor/2014-February/057887.html.  The method of modeling is a random effects approach in which the intra-donor correlation is incorporated into the covariance matrix instead of the linear predictor.   And though as Gordon Smyth states both are good method and the twoway anova approach makes fewer assumptions, the random effects approach is statistically more powerful.  
#
# We have a balanced design in which all donors receive all stimuli (which is really in healthy human donors, life and all of its factors!) Our measurement has so many points -- we are measuring in the skipped exon approach, 42,611 junctions!   It is not possible to encorporate those measurements into the linear predictor.  A two-way ANOVA approach is virtually as powerful as the random effects approach 
# and hence is preferable as it makes fewer assumptions.
#
# For an unbalanced design in which each donor receives only a subset of the stimula, the random effects approach is more powerful.
#
# Random effects approach is equivalent to The first method is twoway anova, a generalization of a paired analysis.
#
# The ijc and sjc are expressionSet objects:
#  counts are obtained from exprs(ijc) and exprs(sjc)
#  
#  operations at the main object level will ripple to the phenotype and expression set information.
#  so filtering occurs on the expressionSet object ijc and sjc

# ### 2.4 Model execution

execute_model <- function (plot, dup, tissue_of_interest, splice_type, fromGTF, tissue_list, ijc, sjc) {

    pdf_sub_directory = '../pdf/'
    csv_sub_directory = '../data/'
    tissue_true   <- tissue_list == tissue_of_interest

    message("Limiting phenotype data to tissue of interest: ",tissue_of_interest)    
    keep = pData(ijc)$SMTSD== tissue_of_interest
       
    ijc <- ijc[,keep==TRUE]
    sjc <- sjc[,keep==TRUE]

    message("dimensions of the pData ijc_tissue: ",dim(pData(ijc)))
    message("dimensions of the pData sjc_tissue: ", dim(pData(sjc)))
    message("dimensions of the exprs ijc_tissue : ", dim(exprs(ijc)))
    message("dimensions of the exprs sjc_tissue: ", dim(exprs(sjc)))

    #
    #  We need to filter low CPM.  This is done ins a similar manner as is done with
    #  genes in an RNA-seq experements.  But rather than ignoring genes, we will ignore
    #  the junctions that do not show significant
    #  levels of expression, this filtering is dependent on minimum CPM: This is the counts per 
    #  million that a gene must have in at least some specified number of samples.
    #        
    #   
    #  number of junctions is given in the first dimension of exprs(ijc) matrix
    #
    #
    #  we keep only junctions for each sex, that have either a minimum number of counts in
    #  both the sjc or ijc)
    #
    count_threshold = min(table(pData(ijc)$SEX))*0.25
    
    male.skip.rs = rowSums(cpm(exprs(sjc)[,pData(sjc)$SEX==1])>=1)
    
    male.inc.rs = rowSums(cpm(exprs(ijc)[,pData(ijc)$SEX==1])>=1)
    
    female.skip.rs = rowSums(cpm(exprs(sjc)[,pData(sjc)$SEX==2])>=1)
    
    female.inc.rs = rowSums(cpm(exprs(ijc)[,pData(ijc)$SEX==2])>=1)
    
    keep= (male.skip.rs>=count_threshold)  & (male.inc.rs>=count_thershold) & 
    
    (female.skip.rs>=count_threshold) & (female.inc.rs>=count_threshold)
    

    
    ijc <- ijc[keep==TRUE,]
    sjc <- sjc[keep==TRUE,]
    
    message("dim pData(ijc)\n",
           paste(dim(pData(ijc))), collapse = " ")
    message("dim exprs(ijc)\n",
           paste(dim(exprs(ijc))), collapse = " ")
    message("dim pData(sjc)\n",
           paste(dim(pData(sjc))), collapse = " ")
    message("dim exprs(sjc)\n",
           paste(dim(exprs(sjc))), collapse = " ")

    
    message("Executing model \n")
    sample_names <- as.character(pData(ijc)$SAMPID)
    # we will add donor as a blocking parameter
    # rather than sample name -- we should use donor for real
    sample     <- factor(sample_names)

    donor    <- rep(sample, 2)
    message("\ndonor size", length(donor)) 

    as_matrix <- cbind(exprs(ijc),exprs(sjc))
    message("dim as_matrix", dim(as_matrix))
            
    message("sex samples: ", table(pData(ijc)$SEX))
    sex2      <- factor(c(rep(pData(ijc)$SEX,2)))
    message("sex samples:\n",
        paste(table(sex2)), collapse = " ")
    message("length sex2: ", length(sex2)) 

    as_event  <- c(rep("ijc",dim(ijc)[2]), rep("sjc", dim(sjc)[2]))
    as_event  <- factor(as_event, levels=c("ijc", "sjc"))
    message("length as_event: ", length(as_event))

    design    <- model.matrix( ~ sex2 + as_event + sex2*as_event)

    colnames(design) <- c("intercept","sex", "as_event","sex*as_event")
    message("\ndim design <- model.matrix( ~sex + as_event + sex*as_event)\n", 
        paste(head(design), collapse = "\n") )

    y <- DGEList(counts=as_matrix, group = sex2)
    y <- calcNormFactors(y, method="RLE")
    y_voom <- voom (y, design=design, plot = plot)

    if (dup==TRUE) {
        dup_cor <- duplicateCorrelation(y_voom$E, design=design, ndups=2, block=donor, weights=y$samples$norm.factors)
        dup_cor$consensus.correlation 
        y_dup_voom <- voom (y, design=design, plot = plot, block = donor, correlation = dup_cor$consensus.correlation) 
    }
    
    Gender <- ifelse(pData(ijc)$SEX==1,"m","f")
    message("\nGenders new size\n", 
        paste(length(Gender), collapse = " ") )
    message("\nplotting y for ijc portion of design <- model.matrix( ~sex + as_event + sex*as_event\n")
     # print the combined exploratory plot
    filename <- paste0(paste0(paste0(pdf_sub_directory, splice_type),
                              snakecase::to_snake_case(tissue_of_interest)),"-y-ijc-MDSplot-100.pdf")
    pdf (filename)
        par(cex=1.5)
        plotMDS(y[,c(1:dim(ijc)[2])], labels=Gender, top=100, col=ifelse(Gender=="m","blue","red"), 
            gene.selection="common")
    dev.off()
    message("\nplotting y_voom for ijc portion of design <- model.matrix( ~sex + as_event + sex*as_event\n")
    filename <- paste0(paste0(paste0(pdf_sub_directory, splice_type),
                              snakecase::to_snake_case(tissue_of_interest)),"-y-voom-ijc-MDSplot-100.pdf")
    pdf (filename)
        par(cex=1.5)
        plotMDS(y_voom[,c(1:dim(ijc)[2])], labels=Gender, top=100, col=ifelse(Gender=="m","blue","red"), 
            gene.selection="common")
    dev.off()
    if (dup == TRUE) {
        filename <- paste0(paste0(paste0(pdf_sub_directory, splice_type),
                           snakecase::to_snake_case(tissue_of_interest)),"-y-dup-voom-ijc-MDSplot-100.pdf")
        pdf (filename)
            par(cex=1.5)
            plotMDS(y_dup_voom[,c(1:dim(ijc)[2])], labels=Gender, top=100, col=ifelse(Gender=="m","blue","red"), 
                gene.selection="common")
        dev.off()
    }
    message("\nplotting y for sjc portion of design <- model.matrix( ~sex + as_event + sex*as_event\n")
    filename <- paste0(paste0(paste0(pdf_sub_directory, splice_type),
                              snakecase::to_snake_case(tissue_of_interest)),"-y-sjc-MDSplot-100.pdf")
    pdf (filename)
        par(cex=1.5)
        plotMDS(y[,c((dim(ijc)[2]+1)):(dim(ijc)[2]+dim(sjc)[2])], labels=Gender, top=100, col=ifelse(Gender=="m","blue","red"), 
            gene.selection="common")
    dev.off()
    
    if (dup == TRUE) {
        message("\nplotting y_voom for sjc portion of design <- model.matrix( ~sex + as_event + sex*as_event\n")    
        filename <- paste0(paste0(paste0(pdf_sub_directory, splice_type),
                              snakecase::to_snake_case(tissue_of_interest)),"-y-voom-sjc-MDSplot-100.pdf")
        pdf (filename)
            par(cex=1.5)
            plotMDS(y_voom[,c((dim(ijc)[2]+1)):(dim(ijc)[2]+dim(sjc)[2])], labels=Gender, top=100, col=ifelse(Gender=="m","blue","red"), 
                gene.selection="common")
        dev.off()
        filename <- paste0(paste0(paste0(pdf_sub_directory, splice_type),
                              snakecase::to_snake_case(tissue_of_interest)),"-y-dup-voom-sjc-MDSplot-100.pdf")
        pdf (filename)
            par(cex=1.5)
            plotMDS(y_dup_voom[,c((dim(ijc)[2]+1)):(dim(ijc)[2]+dim(sjc)[2])], labels=Gender, top=100, col=ifelse(Gender=="m","blue","red"), 
                gene.selection="common")
        dev.off()
        
        fit <- lmFit(y_dup_voom, design=design, block=donor, correlation = dup_cor$consensus.correlation)
    } else {
        fit <- lmFit(y_voom, design=design)
    }
        
    fit <- eBayes(fit, robust=TRUE)
    
    sex_as_events_results         <- topTable(fit, coef="sex*as_event", number=nrow(y_voom))
    sex_as_events_results_refined <- sex_as_events_results$adj.P.Val <= 0.05 & abs(sex_as_events_results$logFC) >= abs(log2(1.5))

    sex_results                   <- topTable(fit, coef="sex", number=nrow(y_voom))
    sex_results_refined           <- sex_results$adj.P.Val <= 0.05 & abs(sex_results$logFC) >= abs(log2(1.5))

    sex_as_events_rnResults <- rownames(sex_as_events_results)
    sex_rnResults           <- rownames(sex_results)
    head(sex_as_events_rnResults)
    head(sex_rnResults)
    head(fromGTF[sex_as_events_rnResults,])

    # use the junctionIDs to get the annotations
    sex_as_events_resultsAnnotations      <- fromGTF[sex_as_events_rnResults,]
    sex_resultsAnnotations                <- fromGTF[sex_rnResults,]
    head(sex_as_events_resultsAnnotations)
    head(sex_resultsAnnotations)
    
    sex_as_events_results_refinedAnnotations<- sex_as_events_resultsAnnotations[sex_as_events_results_refined==TRUE,]
    sex_results_refinedAnnotations          <- sex_resultsAnnotations          [sex_results_refined          ==TRUE,]
    head(sex_as_events_results_refinedAnnotations)
    head(sex_results_refinedAnnotations)

    # geneSymbols are in the annotations 
    sex_as_events_geneSymbols         <- sex_as_events_resultsAnnotations$geneSymbol
    sex_as_events_refined_geneSymbols <- sex_as_events_results_refinedAnnotations$geneSymbol
    sex_geneSymbols                   <- sex_resultsAnnotations$geneSymbol
    sex_refined_geneSymbols           <- sex_results_refinedAnnotations$geneSymbol


    # adjust the rownames to be the geneSymbols rather than junction IDs
    sex_as_events_results_rn   <- paste(sex_as_events_geneSymbols, sex_as_events_rnResults, sep="-")
    sex_results_rn             <- paste(sex_geneSymbols,           sex_rnResults, sep="-")
    message("\n sex_as_events\n", 
        paste(head(sex_as_events_results_rn), collapse = " ") )
    rownames(sex_as_events_results) <- sex_as_events_results_rn
    rownames(sex_results)           <- sex_results_rn

    sex_as_events_filename         = paste0(paste0(paste0(csv_sub_directory, splice_type),
                                                   snakecase::to_snake_case(tissue_of_interest)),'_AS_model_B_sex_as_events.csv')
    sex_as_events_refined_filename = paste0(paste0(paste0(csv_sub_directory, splice_type),
                                                   snakecase::to_snake_case(tissue_of_interest)),'_AS_model_B_sex_as_events_refined.csv',sep='')
    sex_filename                   = paste0(paste0(paste0(csv_sub_directory, splice_type),
                                                   snakecase::to_snake_case(tissue_of_interest)),'_AS_model_B_sex.csv',sep='')
    sex_refined_filename           = paste0(paste0(paste0(csv_sub_directory, splice_type),
                                                   snakecase::to_snake_case(tissue_of_interest)),'_AS_model_B_sex_refined.csv',sep='')
    sex_as_events_genesFilename    = paste0(paste0(paste0(csv_sub_directory, splice_type),
                                                   snakecase::to_snake_case(tissue_of_interest)),'_AS_model_B_sex_as_events_universe.txt',sep='')
    sex_as_events_refined_genesFilename = paste0(paste0(paste0(csv_sub_directory, splice_type),
                                                   snakecase::to_snake_case(tissue_of_interest)),'_AS_model_B_sex_as_events_gene_set.txt',sep='')
    sex_genesFilename              = paste0(paste0(paste0(csv_sub_directory, splice_type),
                                                   snakecase::to_snake_case(tissue_of_interest)),'_AS_model_B_sex_universe.txt',sep='')
    sex_refined_genesFilename      = paste0(paste0(paste0(csv_sub_directory, splice_type),
                                                   snakecase::to_snake_case(tissue_of_interest)),'_AS_model_B_sex_gene_set.txt',sep='')

    write.table(sex_as_events_results, file = sex_as_events_filename, 
                row.names = T, col.names = T, quote = F, sep = ",")
    write.table(sex_as_events_results[sex_as_events_results_refined,], 
                file = sex_as_events_refined_filename, row.names = T, col.names = T, quote = F, sep = ",")
    write.table(sex_results,           file = sex_filename          , 
                row.names = T, col.names = T, quote = F, sep = ",")
    write.table(sex_results [sex_results_refined          ,], file = sex_refined_filename, 
                row.names = T, col.names = T, quote = F, sep = ",")
    write.table(sex_as_events_geneSymbols, file = sex_as_events_genesFilename, 
                row.names = F, col.names = F, quote = F, sep = ",")
    write.table(sex_as_events_refined_geneSymbols,file = sex_as_events_refined_genesFilename, 
                row.names = F, col.names = F, quote = F, sep = ",")
    write.table(sex_geneSymbols,           file = sex_genesFilename          , 
                row.names = F, col.names = F, quote = F, sep = ",")
    write.table(sex_refined_geneSymbols,          file = sex_refined_genesFilename          , 
                row.names = F, col.names = F, quote = F, sep = ",")

    return(0)
}

# ### 3 Execution of All Tissues and All Splicing Variants
#
# Additional values set to enable this notebook to be executed as a nextflow workflow or to run in place with appropriate settings.
# ### 3.1 parameters Setting
#
# 1. Setting `dup=TRUE` causes lengthy execution times.
#
# 2. Setting `plot=TRUE` can overwhelm the saving capacity within a jupyter-lab notebook - 
#    this sets to print all the voom plots.
#    
# 3. Adjusting `splice_type` will allow you to play with a variety of results
#
#    a. all splice types desired to be run:
#     
#     `splice_list       = c("a3ss_","a5ss_","mxe_","ri_","se_")`
#     
#    b. a subset (leaving out say `splice_type = "se_"` since it is the largest, for example)
#     
#     `splice_list       = c("a3ss_","a5ss_","mxe_","ri_")`

# ### 3.2 MAIN routine

getReleasedGTFAnnotations (destDir <- "../data/")

# +
# 1.2.3 get SRR Accession Metadata (available through dbGaP)
srr_metadata <- getSraRunData (destDir <- "../data/")

gtexPhenoDataObj <- getGTExCorrectedRDS (destDir <- "../data/")

# 1.2.6 get reduced Tissue data
tissue_reduction <- getTissueReduction ( "../assets/tissues.tsv" )

# 1.2.8 get GTEx Phenotype Data accessions SRR from releaase
srr_pdata <- getGTExPhenoDataForSRR ("../data/")

# +
#
# Now the meat of our work
# and to reduce to a subset of entire splice_list
# splice_list <- c("a3ss_","a5ss_","mxe_","ri_","se_")
#
plot         <- TRUE
dup          <- TRUE
splice_list  <- c("a3ss_","a5ss_","mxe_","ri_","se_")
tissue_reduction <- tissue_reduction[tissue_reduction$include==1,]
tissue_list  <- factor(snakecase::to_snake_case(as.character(tissue_reduction$SMTSD)))  
write.table(tissue_list, "../data/tissue_list.csv", sep=",", quote=FALSE, row.names=TRUE, col.names=FALSE)

# parameters that change with each splice type (3)
# 1. fromGTF
# 2. ijc
# 3. sjc
# Could run this as a loop - or rather, using a package [package name]
# run this notebook as a nextflow workflow
# Requirements are that all required input are in a bucket data.tar.gz
# and assets 
# for (tissue_index in 1:length(tissue_list)) {
# -
do_test_run <- function() {
    splice_type = "se_"
    res = splice_list %in% splice_type
    tissue_of_interest  = as.vector(as.character(tissue_list[tissue_index]))
    message ("splice_list does contain ", splice_type, " continuing with processing")
    fromGTF    <- read.table("../data/fromGTF.SE.txt", header=TRUE)
    fromGTF    <- eliminateChrYfromGTF (fromGTF)
    ijc_counts <- makeCountsMatrix( filename_gz  <- "../data/rmats_final.se.jc.ijc.txt.gz")
    ijc        <- makeSplicingExpressionSetObject (srr_pdata <- srr_pdata, counts <- ijc_counts)
    ijc        <- reduceSampleSet(tissue_reduction <- tissue_reduction, es <- ijc)
    ijc        <- eliminateChrYfromExpressionSet (fromGTF <- fromGTF, es <- ijc)
    sjc_counts <- makeCountsMatrix( filename_gz  <- "../data/rmats_final.se.jc.sjc.txt.gz")
    sjc        <- makeSplicingExpressionSetObject (srr_pdata <- srr_pdata, counts <- sjc_counts)
    sjc        <- reduceSampleSet(tissue_reduction <- tissue_reduction, es <- sjc)
    sjc        <- eliminateChrYfromExpressionSet (fromGTF <- fromGTF, es <- sjc)
    execute_model (plot, 
                    dup, 
                    tissue_of_interest, 
                    splice_type, 
                    fromGTF, 
                    tissue_list, 
                    ijc, 
                    sjc)
    }
do_test_run()
message("finished test run")     

for (tissue_index in (tissue_list)) {
    # a3ss
    splice_type = "a3ss_"
    res = splice_list %in% splice_type
    tissue_of_interest  = as.vector(as.character(tissue_list[tissue_index]))

    if (sum(res) == 1) {
        message ("splice_list does contain\n",
             paste(splice_type), " continuing with processing\n")
        fromGTF    <- read.table("../data/fromGTF.A3SS.txt", header=TRUE)
        fromGTF    <- eliminateChrYfromGTF (fromGTF)
        ijc_counts <- makeCountsMatrix( filename_gz  <- "../data/rmats_final.a3ss.jc.ijc.txt.gz")
        ijc        <- makeSplicingExpressionSetObject (srr_pdata <- srr_pdata, counts <- ijc_counts)
        ijc        <- reduceSampleSet(tissue_reduction <- tissue_reduction, es <- ijc)
        ijc        <- eliminateChrYfromExpressionSet (fromGTF <- fromGTF, es <- ijc)
        sjc_counts <- makeCountsMatrix( filename_gz  <- "../data/rmats_final.a3ss.jc.sjc.txt.gz")
        sjc        <- makeSplicingExpressionSetObject (srr_pdata <- srr_pdata, counts <- sjc_counts)
        sjc        <- reduceSampleSet(tissue_reduction <- tissue_reduction, es <- sjc)
        sjc        <- eliminateChrYfromExpressionSet (fromGTF <- fromGTF, es <- sjc)
        execute_model (plot, 
                                 dup, 
                                 tissue_of_interest, 
                                 splice_type, 
                                 fromGTF, 
                                 tissue_list, 
                                 ijc, 
                                 sjc)
 
    }
    # a5ss
    splice_type = "a5ss_"
    res = splice_list %in% splice_type
    tissue_of_interest  = as.vector(as.character(tissue_list[tissue_index]))

    if (sum(res) == 1) {
        message ("splice_list does contain\n",
             paste(splice_type), " continuing with processing\n")
        fromGTF    <- read.table("../data/fromGTF.A5SS.txt", header=TRUE)
        fromGTF    <- eliminateChrYfromGTF (fromGTF)
        ijc_counts <- makeCountsMatrix( filename_gz  <- "../data/rmats_final.a5ss.jc.ijc.txt.gz")
        ijc        <- makeSplicingExpressionSetObject (srr_pdata <- srr_pdata, counts <- ijc_counts)
        ijc        <- reduceSampleSet(tissue_reduction <- tissue_reduction, es <- ijc)
        ijc        <- eliminateChrYfromExpressionSet (fromGTF <- fromGTF, es <- ijc)
        sjc_counts <- makeCountsMatrix( filename_gz  <- "../data/rmats_final.a5ss.jc.sjc.txt.gz")
        sjc        <- makeSplicingExpressionSetObject (srr_pdata <- srr_pdata, counts <- sjc_counts)
        sjc        <- reduceSampleSet(tissue_reduction <- tissue_reduction, es <- sjc)
        sjc        <- eliminateChrYfromExpressionSet (fromGTF <- fromGTF, es <- sjc)
        execute_model  (plot, 
                                 dup, 
                                 tissue_of_interest, 
                                 splice_type, 
                                 fromGTF, 
                                 tissue_list, 
                                 ijc, 
                                 sjc)
 
    }
    # mxe
    splice_type = "mxe_"
    res = splice_list %in% splice_type
    tissue_of_interest  = as.vector(as.character(tissue_list[tissue_index]))

    if (sum(res) == 1) {
        message ("splice_list does contain\n",
             paste(splice_type), " continuing with processing\n")
        fromGTF    <- read.table("../data/fromGTF.MXE.txt", header=TRUE)
        fromGTF    <- eliminateChrYfromGTF (fromGTF)
        ijc_counts <- makeCountsMatrix( filename_gz  <- "../data/rmats_final.mxe.jc.ijc.txt.gz")
        ijc        <- makeSplicingExpressionSetObject (srr_pdata <- srr_pdata, counts <- ijc_counts)
        ijc        <- reduceSampleSet(tissue_reduction <- tissue_reduction, es <- ijc)
        ijc        <- eliminateChrYfromExpressionSet (fromGTF <- fromGTF, es <- ijc)
        sjc_counts <- makeCountsMatrix( filename_gz  <- "../data/rmats_final.mxe.jc.sjc.txt.gz")
        sjc        <- makeSplicingExpressionSetObject (srr_pdata <- srr_pdata, counts <- sjc_counts)
        sjc        <- reduceSampleSet(tissue_reduction <- tissue_reduction, es <- sjc)
        sjc        <- eliminateChrYfromExpressionSet (fromGTF <- fromGTF, es <- sjc)
        execute_model  (plot, 
                                 dup, 
                                 tissue_of_interest, 
                                 splice_type, 
                                 fromGTF, 
                                 tissue_list, 
                                 ijc, 
                                 sjc)
 
    }
    # ri
    splice_type = "ri_"
    res = splice_list %in% splice_type
    tissue_of_interest  = as.vector(as.character(tissue_list[tissue_index]))

    if (sum(res) == 1) {
        message ("splice_list does contain\n",
             paste(splice_type), " continuing with processing\n")
        fromGTF    <- read.table("../data/fromGTF.RI.txt", header=TRUE)
        fromGTF    <- eliminateChrYfromGTF (fromGTF)
        ijc_counts <- makeCountsMatrix( filename_gz  <- "../data/rmats_final.ri.jc.ijc.txt.gz")
        ijc        <- makeSplicingExpressionSetObject (srr_pdata <- srr_pdata, counts <- ijc_counts)
        ijc        <- reduceSampleSet(tissue_reduction <- tissue_reduction, es <- ijc)
        ijc        <- eliminateChrYfromExpressionSet (fromGTF <- fromGTF, es <- ijc)
        sjc_counts <- makeCountsMatrix( filename_gz  <- "../data/rmats_final.ri.jc.sjc.txt.gz")
        sjc        <- makeSplicingExpressionSetObject (srr_pdata <- srr_pdata, counts <- sjc_counts)
        sjc        <- reduceSampleSet(tissue_reduction <- tissue_reduction, es <- sjc)
        sjc        <- eliminateChrYfromExpressionSet (fromGTF <- fromGTF, es <- sjc)
        execute_model  (plot, 
                                 dup, 
                                 tissue_of_interest, 
                                 splice_type, 
                                 fromGTF, 
                                 tissue_list, 
                                 ijc, 
                                 sjc)
 
    }
    # se
    splice_type = "se_"
    res = splice_list %in% splice_type
    tissue_of_interest  = as.vector(as.character(tissue_list[tissue_index]))

    if (sum(res) == 1) {
        message ("splice_list does contain ", splice_type, " continuing with processing")
        fromGTF    <- read.table("../data/fromGTF.SE.txt", header=TRUE)
        fromGTF    <- eliminateChrYfromGTF (fromGTF)
        ijc_counts <- makeCountsMatrix( filename_gz  <- "../data/rmats_final.se.jc.ijc.txt.gz")
        ijc        <- makeSplicingExpressionSetObject (srr_pdata <- srr_pdata, counts <- ijc_counts)
        ijc        <- reduceSampleSet(tissue_reduction <- tissue_reduction, es <- ijc)
        ijc        <- eliminateChrYfromExpressionSet (fromGTF <- fromGTF, es <- ijc)
        sjc_counts <- makeCountsMatrix( filename_gz  <- "../data/rmats_final.se.jc.sjc.txt.gz")
        sjc        <- makeSplicingExpressionSetObject (srr_pdata <- srr_pdata, counts <- sjc_counts)
        sjc        <- reduceSampleSet(tissue_reduction <- tissue_reduction, es <- sjc)
        sjc        <- eliminateChrYfromExpressionSet (fromGTF <- fromGTF, es <- sjc)
        execute_model  (plot, 
                                 dup, 
                                 tissue_of_interest, 
                                 splice_type, 
                                 fromGTF, 
                                 tissue_list, 
                                 ijc, 
                                 sjc)
 
    } # for se
}  # for all tissues


# ### Appendix - Metadata
#
# For replicability and reproducibility purposes, we also print the following metadata:
#
# 1. Checksums of **'artefacts'**, files generated during the analysis and stored in the folder directory **`data`**
# 2. List of environment metadata, dependencies, versions of libraries using `utils::sessionInfo()` and [`devtools::session_info()`](https://devtools.r-lib.org/reference/session_info.html)

# ### Appendix - 1. Checksums with the sha256 algorithm

# +
rm (notebookid)
notebookid   = "differentialSplicingJunctionExpressionAnalysis"
notebookid

message("Generating sha256 checksums of the artefacts in the `..data/` directory .. ")
system(paste0("cd ../data && find . -type f -exec sha256sum {} \\;  >  ../metadata/", notebookid, "_sha256sums.txt"), intern = TRUE)
message("Done!\n")

paste0("../metadata/", notebookid, "_sha256sums.txt")

data.table::fread(paste0("../metadata/", notebookid, "_sha256sums.txt"), header = FALSE, col.names = c("sha256sum", "file"))
# -

# ### Appendix - 2. Libraries Metadata

# +
dev_session_info   <- devtools::session_info()
utils_session_info <- utils::sessionInfo()

message("Saving `devtools::session_info()` objects in ../metadata/devtools_session_info.rds  ..")
saveRDS(dev_session_info, file = paste0("../metadata/", notebookid, "_devtools_session_info.rds"))
message("Done!\n")

message("Saving `utils::sessionInfo()` objects in ../metadata/utils_session_info.rds  ..")
saveRDS(utils_session_info, file = paste0("../metadata/", notebookid ,"_utils_info.rds"))
message("Done!\n")

dev_session_info$platform
dev_session_info$packages[dev_session_info$packages$attached==TRUE, ]
# -


