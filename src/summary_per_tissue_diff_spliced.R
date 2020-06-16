# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.1
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# # Supplementary: Counts of differentially spliced exons per examined tissue
#
# This notebook aggregates the results from the differential splicing from (**see** [figureXXX.ipynb](figure1.ipynb)), and more specifically the `limma::topTable()` output dataframes across all tissues in the GTEX cohort and generates summary statistics for the number of genes found to be statistically up or downregulated between male and female subjects.

#  ---
#  
#  **Running this notebook**:
#  
# A few steps are needed before you can run this document on your own. The GitHub repository (https://github.com/TheJacksonLaboratory/sbas) of the project contains detailed instructions for setting up the environment in the **`dependencies/README.md`** document. Before starting with the analysis, make sure you have first completed the dependencies set up by following the instructions described there. If you have not done this already, you will need to close and restart this notebook before running it.
#
# All paths defined in this Notebook are relative to the parent directory (repository). 
#
#  ---
#

# # Loading dependencies

library(dplyr)
library(tidyr)
library(reshape)
library(ggplot2)
library(piggyback)
library(snakecase)

# # Retrieving the results from the Differential Splicing Events using [`ropensci/piggyback`](https://github.com/ropensci/piggyback)
#
# This notebook requires as input data the limma `topTable()` objects from the Differential Gene Expression analysis (see [figureXXX.ipynb]()). We have archived the results from the notebook that generates the results using the method described by the author of the R package [`ropensci/piggyback`](https://github.com/ropensci/piggyback). We use the release named `diff_spliced` (Differential Spliced Exons) in the repo and can be accessed at [TheJacksonLaboratory/sbas/releases/tag/diff_spliced](https://github.com/TheJacksonLaboratory/sbas/releases/tag/diff_spliced). 
#
# For using the [`ropensci/piggyback`](https://github.com/ropensci/piggyback) with private repositories, it is required that a `GITHUB_TOKEN` is stored as a variable in the R environment in which one is working. To generate such a token with sensible default permissions, the R package [usethis]() has a convenient function 
#
# ```R
# # intall.packages("usethis")
# usethis::browse_github_token()
# ```
#
# This will redirect you to GitHub to create your own GitHub token. Once you have the token, you can use it to set up `.Renviron` by typing the following:
#
# ```R
# Sys.setenv(GITHUB_TOKEN = "youractualtokenindoublequotes")
# ```
#
# Then you sre ready to use the function [`piggyback::pb_download()`](https://docs.ropensci.org/piggyback/reference/pb_download.html) to retrieve the `dge.tar.gz` that contains the topTable objects written as .csv file for all 46 examined GTEX tissue cohorts.
#
# ---
#
# ***NOTE***
#
# Avoid using the `.token` argument to share your token directly in the function as you might forget and push your code, along with your private GITHUB_TOKEN to GitHub. If that happens by mistake, it is advised you invalidate the token that has been exposed by accessing [this link](https://github.com/settings/tokens) and clicking `Delete`.
#
# ---

#usethis::browse_github_token()
# this does not seem to work -- did a conda install -c r r-piggyback in a terminal window
#devtools::install_github("ropensci/piggyback@87f71e8", upgrade="never")
Sys.setenv(GITHUB_TOKEN = "ca4182dec71bfe5ff2128e18f4a66a08c0990e71")

# +
piggyback::pb_download(file = "results_sex_as_events.tar.gz",
                       dest = "../data",
                       repo = "TheJacksonLaboratory/sbas",
                       tag  = "rMATs.3.2.5.gencode.V30.GTEx.V8.DE.v.1.0",
                       show_progress = TRUE)

system("mv ../data/results_sex_as_events.tar.gz  ../data/diff_spliced.tar.gz", intern = TRUE)

tar_gz_top_table_diff_splicing_archive <- ("../data/diff_spliced.tar.gz")
if (!file.exists( tar_gz_top_table_diff_splicing_archive )) {
    
    message( paste0("Fetching ", basename(tar_gz_top_table_diff_splicing_archive) , " from GitHub .."))
    # Download archive from GitHub release with tag "dge"
    piggyback::pb_download(file = "diff_spliced.tar.gz",
                           dest = "../data",
                           repo = "TheJacksonLaboratory/sbas",
                           tag  = "diff_spliced",
                           show_progress = TRUE)
    message("Done!\n")
    
    message("Decompressing archive into folder ../data/diff_spliced ..")
    # Decompress in a folder tmp named dge
    system("mkdir -p ../data/diff_spliced && tar xvzf ../data/diff_spliced.tar.gz -C ../data/diff_spliced/", intern = TRUE)
    message("Done!\n")
}
if (file.exists("../data/diff_spliced.tar.gz ")) {
    message("File diff_spliced.tar.gz already available in ../data/ !\n")
    message("Decompressing archive into folder ../data/diff_spliced ..")
    # Decompress in a folder tmp named dge
    system("mkdir -p ../data/diff_spliced && tar xvzf ../data/diff_spliced.tar.gz -C ../data/diff_spliced/", intern = TRUE)
    message("Done!\n")
}





# +
#https://github.com/TheJacksonLaboratory/sbas/releases/download/rMATs.3.2.5.gencode.V30.GTEx.V8.DE.v.1.0/results_dge_ijc_sex.tar.gz

if (file.exists("../data/results_dge_ijc_sex.tar.gz ")) {
    message("File results_dge_ijc_sex.tar.gz already available in ../data/ !\n")
    message("Decompressing archive into folder ../data/results_dge_ijc_sex ..")
    # Decompress in a folder tmp named dge
    system("tar xvzf ../data/results_dge_ijc_sex.tar.gz", intern = TRUE)
    message("Done!\n")
} else {
    message( paste0("Fetching results_dge_ijc_sex.tar.gz", "from GitHub .."))
    # Download archive from GitHub release with tag "dge"
    piggyback::pb_download(file = "results_dge_ijc_sex.tar.gz",
                           dest = "../data",
                           repo = "TheJacksonLaboratory/sbas",
                           tag  = "rMATs.3.2.5.gencode.V30.GTEx.V8.DE.v.1.0",
                           show_progress = TRUE)
    message("Done!\n")
    message("Decompressing archive into folder ../data/results_dge_ijc_sex ..")
    # Decompress in a folder tmp named dge
    system("tar xvzf ../data/results_dge_ijc_sex.tar.gz", intern = TRUE)
    message("Done!\n")
}    
 



# +

# pData_yarn_downloadGTExV8.csv
if (! (file.exists("../data/pData_yarn_downloadGTExV8.csv"))) {
    system("mkdir -p ../data", intern = TRUE)
    message("Fetching pData_yarn_downloadGTExV8.csv from GitHub ..")
    # Download archive from GitHub release with tag "dge"
    piggyback::pb_download(file = "pData_yarn_downloadGTExV8.csv",
                           dest = "../data",
                           repo = "TheJacksonLaboratory/sbas",
                           tag  = "GTExV8.v1.0",
                           show_progress = TRUE)
    message("Done!\n")
}

# fromGTF.tar.gz
if (! (file.exists("../data/fromGTF.tar.gz"))) {
    system("mkdir -p ../data", intern = TRUE)
    message("Fetching fromGTF.tar.gz from GitHub ..")
    # Download archive from GitHub release with tag "dge"
    piggyback::pb_download(file = "fromGTF.tar.gz",
                           dest = "../data",
                           repo = "adeslatt/sbas_gtf",
                           tag  = "rMATS.3.2.5.gencode.v30",
                           show_progress = TRUE)
    message("Done!\n")
    message("Decompressing fromGTF.tar.gz into ../data")
    system("mkdir -p ../data && tar xvfz ../data/fromGTF.tar.gz -C ../data", intern = TRUE)
    message("Done!\n")
    message("Decompressing fromGTF.*.txt.gz into ../data")
    system("gunzip  ../data/fromGTF*.txt.gz ", intern = TRUE)
    message("Done!\n")
}

# -

# # Create a list of named dataframes with the Differentially Spliced Exons `limma::topTable()`s
#
# We will iterate over the list of named dataframes to collect summary statistics. More specifically, retrieve the count of:
# - upregulated
# - downregulated
# - non significant
#
# differentially spliced exons for the contrast males-females per tissue.

# +
suffix_pattern   <- "*_DGE_ijc_sex.csv"
tables_folder    <- "../data/results_dge_ijc_sex/"

tables_filepaths <- list.files(tables_folder, pattern = suffix_pattern, full.names = TRUE)
tables_filenames <- list.files(tables_folder, pattern = suffix_pattern, full.names = FALSE)
# -

all_topTables <- lapply(tables_filepaths,read.csv)
names(all_topTables) <- gsub("_DGE_ijc_sex.csv","", tables_filenames, fixed = TRUE)

# The list named `all_topTables` is the object that holds all the topTable dataframes from each tissue comparison:

length(all_topTables)

summary(all_topTables)

head(all_topTables[[1]] , 21)

# # Loading the feature annotation tab seperated file `fromGTF.SE.txt`

fromGTF.SE <- read.csv("../data/fromGTF.SE.txt",sep = "\t")
fromGTF.SE[["ID"]] <- as.character(fromGTF.SE[["ID"]] )
head(fromGTF.SE)

# # Example with one topTable before iterating over all tissues

# Example topTable and name
topTable <- all_topTables[[21]]
name     <- names( all_topTables)[21]
name


# ## Defining the thresholds for the double criterion filtering:
#
# Criteria:
# - Adjusted p-value < `p_value_cuttoff`
# - Absolute FoldChange > `absFold_change_threshold`

# ----
#
# ***NOTE***
#
# Defining higher in males or females based on the limma design matrix.
# As we have used 1 for encoding the females and 2 for the males, our *reference level* for the contrast in the expression between males and females is 1, the females.
#
#
# From the `limma` documentation:
# >The level which is chosen for the *reference level* is the level which is contrasted against. By default, this is simply the first level alphabetically. We can specify that we want group 2 to be the reference level by either using the relevel function [..]
#
# By convention, we could say that genes with positive log fold change, are higher in males, whereas the opposite holds true for the ones that are observed to have negative log fold change. 
#
# ---

adj.P.Val_threshold  <- 0.05
absFoldChange_cutoff <- 1.5
absFoldChange_cutoff

# Replacing potential `NA` values in the `P.Value`, `adj.P.Val` to keep the columns numeric and avoid coersion.

# replacing NA p-values with p-value = 1
topTable$P.Value[is.na(topTable$P.Value)]     <- 1; 
topTable$adj.P.Val[is.na(topTable$adj.P.Val)] <- 1;

# +
# Add helper variable dummy `FoldChange` variable. Use 2 as base of log, because this is the default from limma
# The following statement calculates a dummy fold change (how many times higher or lower)
# The minus symbol is a convention symbol only! to express eg. a fold change of 0.25 as -4, 4 times lower
topTable$FoldChange_dummy    <-   ifelse(topTable$logFC > 0, 2 ^ topTable$logFC, -1 / (2 ^ topTable$logFC))                    

# Add helper variable `abs_logFC`.
topTable$abs_logFC <- abs(topTable$logFC)

# Add helper variable `abundance` for up, down, non_signif
topTable$abundance                                                  <- "non_signif"
topTable$abundance[ ((topTable$logFC >=   log2(absFoldChange_cutoff)) & (topTable$adj.P.Val <= adj.P.Val_threshold )) ]   <- "higher"
topTable$abundance[ ((topTable$logFC <=  -log2(absFoldChange_cutoff)) & (topTable$adj.P.Val <= adj.P.Val_threshold )) ]   <- "lower"
# -



dim(topTable[ topTable$abundance == "non_signif" ,])
dim(topTable[ topTable$abundance == "higher" ,])
dim(topTable[ topTable$abundance == "lower" ,])

table(topTable$abundance)

# # Define a vector with the columns to keep in the annotated from GTF `topTable` object

#toKeep <- c("Geneid","logFC","FoldChange_dummy", "adj.P.Val", "abundance")
toKeep <- colnames(topTable)

head(topTable[ , colnames(topTable) %in% toKeep ],2)

# # Annotate topTable with the fromGTF file

topTable$junctionGeneSymbol <- rownames(topTable)
topTable$symbol  <- sub("-[^-]+$","", topTable$junctionGeneSymbol)
topTable$ID <- sub("^.+-", "", topTable$junctionGeneSymbol)

head(topTable, 2)

name
dim(topTable)
dim(topTable  [ topTable$abundance != "non_signif",  ])
dim(topTable  [ ((abs(topTable$logFC > log2(absFoldChange_cutoff)) )  & (topTable$adj.P.Val <= adj.P.Val_threshold )) ,  ])
dim(topTable  [ ((abs(topTable$logFC < -log2(absFoldChange_cutoff)) )  & (topTable$adj.P.Val <= adj.P.Val_threshold )) ,  ])

expression_abundance <- t(table(topTable$abundance))
expression_abundance

expression_abundance <- t(table(topTable$abundance))
signif <- as.data.frame.matrix(expression_abundance)

signif

annotatedTopTable <- dplyr::left_join(topTable, fromGTF.SE, by = "ID")

table(annotatedTopTable$abundance)

head(annotatedTopTable)

# To avoid errors in the cases that we might have none lower or none higher, and the matrix might be missing columns we will create a template data.frame and also add the column that might be missing if lower or higher genes is equal to 0.

signif_template <- structure(list(higher = integer(0), 
                                   lower = integer(0), 
                                   non_signif = integer(0)), 
                              row.names = integer(0), class = "data.frame")
signif_template

# In the for-loop we will check if both columns `lower`, `higher` are present, if not add the column and zero count to create the expected shape of the dataframe:
#
# ```R
# signif <- as.data.frame.matrix(expression_abundance)
# if(! ("higher" %in% colnames(signif))) { 
#     
#     signif$higher <- 0
# }
# if(! ("lower" %in% colnames(signif))) { 
#
#     signif$lower <- 0
# }
# ```

# Now we can add some more summary statistics eg percentage of genes lower, higher or non-significantly different, 

signif$tissue <- name
signif$sum    <- signif$non_signif + signif$higher + signif$lower
toKeepInOrder <- c("tissue", "non_signif", "lower", "higher", "% lower", "% higher", "% non-signif")
signif$`% higher`     <-  round(signif$higher / signif$sum  * 100, 2)
signif$`% lower`      <-  round(signif$lower / signif$sum  * 100, 2)
signif$`% non-signif` <-  round(signif$non_signif / signif$sum  * 100, 2)
signif <- signif[, toKeepInOrder]
signif

# # Summary table of differentially expressed genes between male and female acrosss tissues
#
# Above we demonstrate for one example limma `topTable`. Let's now iterate over all tissue and create an aggregated table of counts of differentially expressed or non-significantly altered between the two sexes.

# +
summary_signif <-structure(list(tissue = character(0), 
                            non_signif = integer(0), 
                            lower = integer(0),
                            higher = integer(0),
                            `% lower` = numeric(0), 
                            `% higher` = numeric(0), 
                            `% non-signif` = numeric(0)), 
                       row.names = integer(0), 
                       class = "data.frame")

signif_template <- structure(list(higher = integer(0), 
                                   lower = integer(0), 
                                   non_signif = integer(0)), 
                              row.names = integer(0), class = "data.frame")

signif_per_tissue <- structure(list(logFC = numeric(0), AveExpr = numeric(0), t = numeric(0), 
                        P.Value = numeric(0), adj.P.Val = numeric(0), B = numeric(0), 
                        initial_gene_id = character(0), gene_id = character(0), abs_logFC = numeric(0), 
                        FoldChange_dummy = numeric(0), abundance = character(0), 
                        GeneSymbol = character(0), Chromosome = character(0), Class = character(0), 
                        Strand = character(0), tissue = character(0)), row.names = integer(0), class = "data.frame")


for (i in seq_along(all_topTables)){
    topTable <- all_topTables[[i]]
    topTable$junctionGeneSymbol <- rownames(topTable)
    topTable$symbol  <- sub("-[^-]+$","", topTable$junctionGeneSymbol)
    topTable$ID <- sub("^.+-", "", topTable$junctionGeneSymbol)  
    topTable$splice_type  <- "SE"
    topTable$featurespace <- "sex_as_events"
    topTable <- dplyr::left_join(topTable, fromGTF.SE, by = "ID")
    name     <- names(all_topTables)[i] 
    # replacing NA p-values with p-value = 1
    topTable$P.Value[is.na(topTable$P.Value)]     <- 1; 
    topTable$adj.P.Val[is.na(topTable$adj.P.Val)] <- 1;
    topTable$abs_logFC <- abs(topTable$logFC)
    # Add helper variable dummy `FoldChange` variable. Use 2 as base of log, because this is the default from limma
    # The following statement calculates a dummy fold change (how many times higher or lower)
    # The minus symbol is a convention symbol only! to express eg. a fold change of 0.25 as -4, 4 times lower
    topTable$FoldChange_dummy    <-   ifelse(topTable$logFC > 0, 2 ^ topTable$logFC, -1 / (2 ^ topTable$logFC))                    

    # Add helper variable `abs_logFC`.
    topTable$abs_logFC <- abs(topTable$logFC)

    # Add helper variable `abundance` for up, down, non_signif
    topTable$abundance                                                  <- "non_signif"
    topTable$abundance[ ((topTable$logFC >=   log2(absFoldChange_cutoff)) & (topTable$adj.P.Val <= adj.P.Val_threshold )) ]   <- "higher"
    topTable$abundance[ ((topTable$logFC <=  -log2(absFoldChange_cutoff)) & (topTable$adj.P.Val <= adj.P.Val_threshold )) ]   <- "lower"
    topTable_signif <- topTable[ topTable$abundance != "non_signif", ]
    topTable_signif$tissue <- name
    signif_per_tissue <- rbind(signif_per_tissue, topTable_signif )
    data.table::fwrite(file = paste0("../data/signif_", snakecase::to_snake_case(name), ".csv"), topTable_signif)
    expression_abundance <- t(table(topTable$abundance))
    signif <- as.data.frame.matrix(expression_abundance)
    if(! ("higher" %in% colnames(signif))) {
        signif$higher <- 0
    }
    if(! ("lower" %in% colnames(signif))) {
        signif$lower <- 0
    }
    signif$tissue <- name
    signif$sum    <-   signif$non_signif + signif$higher + signif$lower
    toKeepInOrder <- c("tissue", "non_signif", "lower", "higher", "% lower", "% higher", "% non-signif")
    signif$`% higher`     <-  round(signif$higher / signif$sum  * 100, 2)
    signif$`% lower`      <-  round(signif$lower / signif$sum  * 100, 2)
    signif$`% non-signif` <-  round(signif$non_signif / signif$sum  * 100, 2)
    signif <- signif[, toKeepInOrder]
    summary_signif <- rbind(summary_signif, signif)   
}
# -

summary_signif <- summary_signif[order(summary_signif$`% non-signif`), ]
head(summary_signif , 2)
head(signif_per_tissue, 2)

# # Defining higher in males or females based on the limma design matrix
# As we have used 1 for encoding the males and 2 for the females, our *reference level* for the contrast in the expression between males and females is 1, the males.
#
#
# From the `limma` documentation:
# >The level which is chosen for the *reference level* is the level which is contrasted against. By default, this is simply the first level alphabetically. We can specify that we want group 2 to be the reference level by either using the relevel function [..]
#
# By convention, we could say that splicing events with positive log fold change, are higher in females, whereas the opposite holds true for the ones that are observed to have negative log folde change. 

summary_signif$`higher in males`   <- summary_signif$lower
summary_signif$`higher in females` <- summary_signif$higher
head(summary_signif[summary_signif$tissue == "artery_aorta", ])

summary_signif$lower

# # Preparing the summary table for plotting
#
# We will need to aggregate the number of genes in one column in order to be able to plot, and also convert the `Tissue` column to a factor. We will use the `reshape` R package to *melt* the dataframe from a wide to a long version, as described above:

toPlot <- summary_signif[, c( "tissue", "higher in males", "higher in females")]
toPlot <- reshape::melt(toPlot, id=c("tissue"))
toPlot$tissue <- as.factor(toPlot$tissue)
colnames(toPlot) <- c("Tissue", "Sex Bias", "Number of Exons")
head(toPlot[toPlot$Tissue == "artery_aorta", ])

head(toPlot)


# +
options(repr.plot.width=8, repr.plot.height=10)

ggplot(toPlot, aes(x = Tissue, y = `Number of Exons`, fill = `Sex Bias`)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_manual (values = c( "higher in males" = "#4A637B" , "higher in females" = "#f35f71")) + 
  theme(text              = element_text(color = "#4A637B", face = "bold", family = 'Helvetica')
        ,plot.caption     = element_text(size =  6, color = "#8d99ae", face = "plain", hjust= 1.05) 
        ,plot.title       = element_text(size =  8, color = "#2b2d42", face = "bold", hjust= 0.5)
        ,axis.text.y      = element_text(angle =  0, size = 6, color = "#8d99ae", face = "bold", hjust=1.1)
        ,axis.text.x      = element_text(angle = 70, size = 6, color = "#8d99ae", face = "bold", hjust=1.1)
        ,axis.title.x     = element_blank()
        ,axis.ticks.x     = element_blank()
        ,axis.ticks.y     = element_blank()
        ,plot.margin      = unit(c(1,1,1,1),"cm")
        ,panel.background = element_blank()
        ,legend.position  = "right") +
  

  geom_text(aes(y = `Number of Exons` + 40, 
                label = `Number of Exons`),
                size = 2,
                color     = "#4A637B",
                position  =  position_dodge(width = 1),
                family    = 'Helvetica') +
  
  labs(title   = "Number of genes with higher expression of included exons in each sex per tissue\n",
       caption = "\nsource: 'The impact of sex on alternative splicing'\n doi: https://doi.org/10.1101/490904",
       y   = "\nNumber of Differentially Expressed Genes")  + coord_flip()


# -

# # Mutually exclusive sex biased genes (higher expression in one or the other sex only)
#
#
# The dataframe `signif_per_tissue` contains all the information for the genes that were significantly higher in either of the two sexes. WLet's examine how many mutually exclusive genes were found across all examined tissues. Ensembl encodes as `Chromosome` the chromosomal position, so we will create the required variables to retrieve only the chromosome information for producing summary statistics.

dput(colnames(signif_per_tissue))
t(head(signif_per_tissue$chr,1))

chromosome_column <- "chr"
signif_per_tissue$Chromosome <-  signif_per_tissue[[chromosome_column]]
signif_per_tissue$higher_in  <- 0
signif_per_tissue$higher_in[(signif_per_tissue$abundance == "lower" )] <- "males"
signif_per_tissue$higher_in[(signif_per_tissue$abundance == "higher" )] <- "females"
toKeepInOrder <- colnames(signif_per_tissue)
signif_per_tissue <- signif_per_tissue[, toKeepInOrder]
head(signif_per_tissue)

# # Examine mutually exclusive genes upregulated in each sex

# +
feature_id_column <- "junctionGeneSymbol"
female_biased <- unique(signif_per_tissue[[feature_id_column]] [ signif_per_tissue$higher_in == "females" ] )
male_biased   <- unique(signif_per_tissue[[feature_id_column]] [ signif_per_tissue$higher_in == "males"  ] )

length(male_biased)
length(female_biased)

# +
## Present in both

length((intersect(male_biased, female_biased)))
length((intersect(female_biased, male_biased)))

intersect <- (intersect(male_biased, female_biased))


# +
## Only in males
length(male_biased[! (male_biased %in% intersect)])

## Only females
length(female_biased[! (female_biased %in% intersect)])

# +
perc_only_male <-  length(male_biased[! (male_biased %in% intersect)]) / length(male_biased) * 100
perc_only_female <-  length(female_biased[! (female_biased %in% intersect)]) / length(female_biased) * 100

head( signif_per_tissue[ signif_per_tissue[[feature_id_column]] %in% male_biased[! (male_biased %in% intersect)],  ] , 4 )

message(round(perc_only_male, 2), " % of differentially spliced exons higher in males only found to be significantly different in males")
message(round(perc_only_female,2), " % of differentially spliced exon higher in females only found to be significantly different in females")
# -

# ## Significantly higher only in males

dim(signif_per_tissue[ signif_per_tissue[[feature_id_column]] %in% male_biased[! (male_biased %in% intersect)],  ])
only_male_genes <- signif_per_tissue[ signif_per_tissue[[feature_id_column]] %in% (male_biased[! (male_biased %in% intersect)]) ,  ]
head(only_male_genes[ order(only_male_genes[[feature_id_column]] ), ], 5)

# +
# See 8.1.1 enquo() and !! - Quote and unquote arguments in https://tidyeval.tidyverse.org/dplyr.html

only_male_genes %>% 
    count( symbol, sort = TRUE) %>%
    head(20)
# -

# ## Significantly higher only in females

# +
only_female_genes <- signif_per_tissue[ signif_per_tissue[[ feature_id_column ]] %in% (female_biased[! (female_biased %in% intersect)]) ,  ]

head(only_female_genes[ order(only_female_genes[[ feature_id_column ]] ), ], 10)
# -

only_female_genes %>% 
    count( symbol, sort = TRUE) %>%
    head(20)

# # Examine number of differentially expressed genes per chromosome per sex

# +
signif_per_tissue$Chromosome <- as.factor(signif_per_tissue$Chromosome)
signif_per_tissue$higher_in <- as.factor(signif_per_tissue$higher_in)

signif_per_tissue %>% 
    group_by(Chromosome,higher_in) %>%  
    count()  -> signif_per_chrom_per_sex
# -

signif_per_chrom_per_sex

# ## Metadata
#
# For replicability and reproducibility purposes, we also print the following metadata:
#
# 1. Checksums of **'artefacts'**, files generated during the analysis and stored in the folder directory **`data`**
# 2. List of environment metadata, dependencies, versions of libraries using `utils::sessionInfo()` and [`devtools::session_info()`](https://devtools.r-lib.org/reference/session_info.html)

# +
notebook_id   = "summary_per_tissue_diff_expressed"

message("Generating sha256 checksums of the artefacts in the `..data/` directory .. ")
system(paste0("cd ../data/ && sha256sum **/*csv > ../metadata/", notebook_id, "_sha256sums.txt"), intern = TRUE)
system(paste0("cd ../data/ && sha256sum *csv >> ../metadata/", notebook_id, "_sha256sums.txt"), intern = TRUE)

message("Done!\n")

data.table::fread(paste0("../metadata/", notebook_id, "_sha256sums.txt"), header = FALSE, col.names = c("sha256sum", "file"))
# -

# ### 2. Libraries metadata

# +
dev_session_info   <- devtools::session_info()
utils_session_info <- utils::sessionInfo()

message("Saving `devtools::session_info()` objects in ../metadata/devtools_session_info.rds  ..")
saveRDS(dev_session_info, file = paste0("../metadata/", notebook_id, "_devtools_session_info.rds"))
message("Done!\n")

message("Saving `utils::sessionInfo()` objects in ../metadata/utils_session_info.rds  ..")
saveRDS(utils_session_info, file = paste0("../metadata/", notebook_id ,"_utils_info.rds"))
message("Done!\n")

dev_session_info$platform
dev_session_info$packages[dev_session_info$packages$attached==TRUE, ]
# -

# # Calculating the sex-biased splicing index
# The normalized sex-biased splicing index is defined as the number of statistically significant splicing events per 1000 exons in the chromosome.

dim(signif_per_chrom_per_sex)

# Sorry, I cannot do this in R. Here is Python (ugly script but works)
#
# import csv
# import gzip
# import re
# from collections import defaultdict
#
# fname = 'Homo_sapiens.GRCh38.100.chr_patch_hapl_scaff.gtf.gz'
#
# chrom2exons = defaultdict(set)
#
# with gzip.open(fname, 'rt') as f:
#     cr = csv.reader(f, delimiter='\t', quotechar='"')
#     for row in cr:
#         #print(row)
#         if row[0].startswith('#'):
#             continue
#         chrom = row[0]
#         annots = row[8]
#         fields = annots.split(";")
#         exon = re.compile(r'exon_id "(ENSE\d+)"')
#         for f in fields:
#             itm = f.strip()
#             match = exon.match(itm)
#             if match:
#                 exonid = match.group(1)
#                 chrom2exons[chrom].add(exonid)
#
# g = open('chrom2exons.txt', 'wt')
# for k, v in chrom2exons.items():
#     print("chr{}: n={}".format(k, len(v)))
#     g.write("{}\t{}\n".format(k, len(v)))
# g.close()
#
# 1	69381
# 2	55599
# 3	46452
# 4	29749
# 5	34789
# 6	33817
# 7	35973
# X	22471
# 8	28489
# 9	26460
# 11	43212
# 10	26514
# 12	42925
# 13	13193
# 14	25994
# 15	28720
# 16	36285
# 17	45142
# 18	13360
# 20	16704
# 19	44166
# Y	2908
# 22	16411
# 21	8830
# MT	37

signif_per_chrom_per_sex

only_female_genes %>% 
    count( !!GENE_ID, GeneSymbol, Class, sort = TRUE) %>%
    head(20)

signif_per_tissue %>% 
     group_by(Chromosome) %>%  
    count()  -> signif_per_chrom

signif_per_chrom

chrom2exon_filename = '../assets/canon_chrom2exons.txt'
if (! file.exists(chrom2exon_filename)) {
    message("Could not find canon_chrom2exons.txt file")
}
c2e_df = read.csv(chrom2exon_filename, sep='\t', header=FALSE)
colnames(c2e_df) <- c("Chromosome","exons")
head(c2e_df) # 25 chromosomes including MT

df2 <- merge(signif_per_chrom, c2e_df, by="Chromosome")
head(df2)

# calculate splicinig index
library(tidyverse)

df2 %>% 
  mutate(Index = 1000 * n/exons) -> df3

df4 <- df3[-25,] # remove the Y chromosome
df4 <- df4[-23,] # remove the MT chromosome

res_sorted <- df4[order(df4$Index, decreasing=TRUE),]
res_sorted


res_sorted$Chromosome <- factor(res_sorted$Chromosome, levels = res_sorted$chr)
res_sorted

# set the colors
npgBlue<- rgb(60/256,84/256,136/256,1)
npgRed <- rgb(220/256,0,0,0.5)
npgGreen <- rgb(0,160/256,135/256,1)
npgBrown <- rgb(126/256,97/256,72/256,1)

# make the plot 
figure2b <- ggplot(res_sorted, aes(x = Chromosome, y = Index, size = n)) +
  geom_point(color=npgBlue) +
  theme_bw() +
  theme(axis.text.x = element_text(size=14, angle = 270, hjust = 0.0, vjust = 0.5),
	axis.text.y = element_text(size=16),
	axis.title.x = element_blank(),
	axis.title.y = element_text(face="plain", colour="black",
                                    size=18),
	legend.title=element_blank(),
	legend.text = element_text(face="plain", colour="black",
                                   size=14)) +
  scale_fill_viridis_c() +
  ylab(paste("Sex-biased splicing index ")) +
  xlab("Chromosomes") +
  guides(size = guide_legend(title = "Number of ASE"))
figure2b



