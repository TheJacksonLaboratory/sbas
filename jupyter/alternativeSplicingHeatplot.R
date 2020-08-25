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

# # **Alternative Splicing  Heatplot**
# This notebook generates a heat plot representing sex-biased differential gene expression as well as a plot showing the counts of differentially expressed genes per tissue.
#
# The values in the heatmap represent the correlation (similarity in the fold-changes) between male and female samples, with the values in the heatmap being the correlation between the vectors of fold changes of the tissues.
#
# The assumptions made before rendering the heatmaps 
# 1. Get differential gene expression (DGE) files
# 2. Use the ``../assets/tissues.tsv`` to limit the tissues to those with at least 50 samples in each sex (``tissues.tsv`` was produced by a Python script)
# 3. Use the pattern for the differentially expressed genes **"../data/*_AS_model_B_sex_as_events.csv"** to get all the values for the matrix.

# ## **Running this notebook**:
#
# See the README for setting up prerequisites for the notebook.

# ## 1. Setup 
#
# Assumes the `countGenesAndEvents.ipynb` notebook was run -- unpacking the results from the differential Gene Expression Analysis as run in the `differentialAlternativeSplicingJunctionAnalysis.ipynb` notebook.

suppressMessages({
    options(warn = -1) 
    library(stringr)
    library(edgeR)
    library(pheatmap)
    library(dplyr)
    library(ggplot2)
    library(scales)
    library(viridis)
    library(magrittr)
    library(ComplexHeatmap)
    library(circlize)
    library(snakecase)
    library(pvclust)
    library(devtools)
    Sys.setenv(TAR = "/bin/tar") # for gzfile
})

# ## 2. Making the matrices

# ### 2.1 Read in all the alternative Splicing Junction Analysis results
#
# Start with the Skipped Exon alternative Splicing Junction Analysis results

results_dir <- "../data/"
pattern     <- "*_AS_model_B_sex_as_events_refined.csv"
files       <- list.files(path = results_dir, pattern = pattern)
message("Number of se AS files found with *AS_model_B_sex_as_events_refined.csv pattern: ",
        length(files))
se_files <- files[grepl("^se_", files)]
length(se_files)

# ### 2.2 read in the curated "../assets/tissues.tsv" 
#
# The **`../assets/tissues.tsv`** file contains an indication to include the tissue if the file has at least **50** samples in that tissue with either **male** or **female** sex reporting

# +
# read in all requirements so that the stage is properly set -- 
# if it is clear here -- it will remain clear for the rest of the time
# tissues.tsv contains the subset of files desired for analysis.
tissue_reduction <- read.table(file="../assets/tissues.tsv", header=TRUE, sep="\t",
                               skipNul=FALSE, stringsAsFactors = FALSE)
colnames(tissue_reduction)  <- c("SMTSD","female","male","include","display_name")
tissue_reduction$SMTSD <- factor(snakecase::to_snake_case(as.character(tissue_reduction$SMTSD)))
# only include those tissues we wish to continue with
table(tissue_reduction$include)
tissue_reduction <- tissue_reduction[tissue_reduction$include==1,]

message("Number of tissues with >=50 samples each in ../assets/tissues.tsv (tissue_reduction)",
        paste(dim(tissue_reduction), collapse=" "))
# -

# ### 2.3 model rownames from fromGTF files
#
# Arbitrarily using the first file, to obtain the ordered rownames for assignment to the remainder of the files in the construction of the matrix.

a3ss_annot <- read.table(file = "../data/fromGTF.A3SS.txt", sep = "\t", quote = "\"", header = T, stringsAsFactors = F)
a5ss_annot <- read.table(file = "../data/fromGTF.A5SS.txt", sep = "\t", quote = "\"", header = T, stringsAsFactors = F)
mxe_annot <- read.table(file = "../data/fromGTF.MXE.txt", sep = "\t", quote = "\"", header = T, stringsAsFactors = F)
ri_annot <- read.table(file = "../data/fromGTF.RI.txt", sep = "\t", quote = "\"", header = T, stringsAsFactors = F)
se_annot <- read.table(file = "../data/fromGTF.SE.txt", sep = "\t", quote = "\"", header = T, stringsAsFactors = F)

# ### 2.4 make a new column combining idx with geneSymbol
#
# all of the model_B toptable events are annotated with **geneSymbol-alternativeSplicingEventID**.  These are alternative Splicing Event specific.   So we will have a matrix per alternative splicing event.

a3ss_annot$geneASID <- as.character(paste(a3ss_annot$geneSymbol, a3ss_annot$ID, sep="-"))
a5ss_annot$geneASID <- as.character(paste(a5ss_annot$geneSymbol, a5ss_annot$ID, sep="-"))
mxe_annot$geneASID  <- as.character(paste(mxe_annot$geneSymbol,  mxe_annot$ID, sep="-"))
ri_annot$geneASID   <- as.character(paste(ri_annot$geneSymbol,   ri_annot$ID, sep="-"))
se_annot$geneASID   <- as.character(paste(se_annot$geneSymbol,   se_annot$ID, sep="-"))


# ### 2.5 Heatmap gene-junction ids are union of events from all tissues
#
# Junctions were removed for statistical requirements -- there are not the same number of results for all of the tissues.  As such we need to build the union of all the ids, order this and use this to enable making all of the matrices.

for (i in 1:length(se_files)) {
    fullfilename  <- paste0("../data/",se_files[i])
    logFC_result  <- read.csv(fullfilename)
    filegenesids  <- as.character(rownames(logFC_result))
    if (i == 1) {
        genesids <- filegenesids
    }
    genesids      <- union(as.character(genesids),filegenesids)
}
length(se_files)
idx <- str_order(genesids)
genesids <- genesids[idx]
length(genesids)
head(genesids)

# ### 2.7 With the union of gene-junction ids make a matrix of logFC and adjusted pValue data

# all these results are already refined, i.e. logFC > log(1.5) and adjusted P-value < 0.05

files <- se_files
logFC_mat <- matrix(0.0, nrow=length(genesids), ncol = length(files))
pVal_mat  <- logFC_mat
    fullfilename  <- paste0("../data/",files[i])
    logFC_result  <- read.csv(fullfilename)
    filegenesids  <- as.character(rownames(logFC_result))

    for (j in 1:length(filegenesids)) {
        genesids_match <- genesids %in% filegenesids[j]
        logFC_mat[genesids_match==TRUE,i] <- logFC_result[filegenesids[j],"logFC"]
        pVal_mat [genesids_match==TRUE,i] <- logFC_result[filegenesids[j],"adj.P.Val"]
    }


# +
files <- se_files
logFC_mat <- matrix(0.0, nrow=length(genesids), ncol = length(files))
pVal_mat  <- logFC_mat
rownames(logFC_mat) <- as.character(genesids)
for (i in 1:length(files)) {
    fullfilename  <- paste0("../data/",files[i])
    logFC_result  <- read.csv(fullfilename)
    filegenesids  <- as.character(rownames(logFC_result))
    # match the filegenesids in the list of unioned genesids for placing the significant gene-junction values
    # in the matrix where rows are gene-junctions and columns are tissues (filegenesids)
    for (j in 1:length(filegenesids)) {
        genesids_match <- genesids %in% filegenesids[j]
        logFC_mat[genesids_match==TRUE,i] <- logFC_result[filegenesids[j],"logFC"]
        pVal_mat [genesids_match==TRUE,i] <- logFC_result[filegenesids[j],"adj.P.Val"]
    }
}

mat_colnames = files
for (i in 1:length(mat_colnames)) {
    mat_colnames[i] <- gsub(pattern,"", mat_colnames[i], fixed = TRUE)
    mat_colnames[i] <- substring(mat_colnames[i],4)
}

display_name = ifelse(mat_colnames == tissue_reduction$SMTSD, 
                                      tissue_reduction$display_name,
                                      mat_colnames)
colnames(logFC_mat) <- as.character(display_name)
colnames(pVal_mat)  <- as.character(display_name)
# -

logFC_mat_rs <- rowSums(abs(logFC_mat))
names(logFC_mat_rs) <- as.character(rownames(logFC_mat))
head(logFC_mat_rs)
plot(log10(logFC_mat_rs))
length(names(logFC_mat_rs))
names(logFC_mat_rs[2678:2683])
logFC_mat[2678:2683,]

# # 2.9 generate a heatmap using agglomeration methods by Ward and distance metric by Minkowski by tissues

colnames(logFC_mat) <- as.character(tissue_reduction$display)

logFC_mat_pv <- pvclust(logFC_mat, method.hclust="ward.D2", method.dist="minkowski", nboot = 1000)

plot(logFC_mat_pv)
pdf("../pdf/logFC_mat_pvclust_ward_minkowski.pdf")
plot(logFC_mat_pv)
dev.off()

# ## the fold change is the qualifying event, but a matrix can be shown of the the junction expression values 
#
# Using the function from the differentialSplicingJunctionAnalysis.ipynb make a counts matrix of the IJC and SJC matrices and show the heatmap of their results

# +
makeCountsMatrix <- function (filename_gz) {
    message("\nloading ", paste(filename_gz), collapse=" ")
    counts <- data.table::fread(filename_gz)
    message("done!")
    rownames(counts) <- counts$ID
    counts <- counts[,-1]
    counts <- data.matrix(counts)
    return(counts)
}

ijc_counts <- makeCountsMatrix( filename_gz  <- "../data/rmats_final.se.jc.ijc.txt.gz")
sjc_counts <- makeCountsMatrix( filename_gz  <- "../data/rmats_final.se.jc.sjc.txt.gz")


# -

logFC_mat[1:10,1:10]


# ### 2.10 Calculate the correlation between the tissues using the logFC
#
# Calculate the correlation distances between the tissues using the logFC within each of the tissues.   Clustering by similar expression.

dist_mat<-as.matrix(cor(logFC_mat))
colnames(dist_mat) <- as.character(tissue_reduction$display)
rownames(dist_mat) <- as.character(tissue_reduction$display)
desc(dist_mat)

sum(is.na(dist_mat))
sum(!is.na(dist_mat))
dist_mat_no_NA <- dist_mat
dist_mat_no_NA[is.na(dist_mat_no_NA)] = 0
sum(is.na(dist_mat_no_NA))
sum(!is.na(dist_mat_no_NA))
sum(dist_mat_no_NA <0)
dist_mat_no_NA_NZ <- dist_mat_no_NA
dist_mat_no_NA_NZ[dist_mat_no_NA_NZ <0] = 0
sum(dist_mat_no_NA_NZ <0)
sum(dist_mat_no_NA_NZ >0)
base_mean = rowMeans(dist_mat_no_NA_NZ)
base_mean

# ## 3. Render the heatmaps of the distance correlations

# ### 3.1 heatmap dist_mat_no_NA_NZ

pheatmap(as.matrix(dist_mat_no_NA_NZ), clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", fontsize = 6)
hm.parameters <- list(dist_mat_no_NA_NZ, fontsize = 6)
do.call("pheatmap", c(hm.parameters,  filename="../pdf/se_alternativeSplicingDistanceCorrelationHeatmapAlllogFC.pdf"))

# ### 3.2 heatmap normalize quantiles dist_mat_no_NA_NZ

dist_mat_NQ <- normalizeQuantiles(dist_mat_no_NA_NZ)
dist_mat <- as.matrix(cor(dist_mat_NQ))
pheatmap(as.matrix(dist_mat), clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", fontsize = 6)
hm.parameters <- list(dist_mat, fontsize = 6)
do.call("pheatmap", c(hm.parameters,  filename="../pdf/se_alternativeSplicingDistanceCorrelationAlllogFC_NQ.pdf"))

# ### 3.3 save the dist_mat rds object

# +
message("Saving dist_mat object")
saveRDS(object = dist_mat, file = "../data/se_as_dist_mat.rds")
message("Done!")
# -

# ## Appendix Metadata
#
# For replicability and reproducibility purposes, we also print the following metadata:
#
# ### Appendix.1. Checksums with the sha256 algorithm
# 1. Checksums of **'artefacts'**, files generated during the analysis and stored in the folder directory **`data`**
# 2. List of environment metadata, dependencies, versions of libraries using `utils::sessionInfo()` and [`devtools::session_info()`](https://devtools.r-lib.org/reference/session_info.html)

figure_id   = "alternativeSplicingHeatmap"

# ### Appendix.2. Libraries

# +
dev_session_info   <- devtools::session_info()
utils_session_info <- utils::sessionInfo()

message("Saving `devtools::session_info()` objects in ../metadata/devtools_session_info.rds  ..")
saveRDS(dev_session_info, file = paste0("../metadata/", figure_id, "_devtools_session_info.rds"))
message("Done!\n")

message("Saving `utils::sessionInfo()` objects in ../metadata/utils_session_info.rds  ..")
saveRDS(utils_session_info, file = paste0("../metadata/", figure_id ,"_utils_info.rds"))
message("Done!\n")

dev_session_info$platform
dev_session_info$packages[dev_session_info$packages$attached==TRUE, ]
# -


