# -*- coding: utf-8 -*-
# # Analysis Notebook - Count Genes and Events
#
# This notebook processes the raw counts as provided by rMATS and performs some descriptive statistical analysis. It is used to produce the following outputs. 
#
# ## Data files created by this notebook
# Output text files are written to the ``data/`` directory (at the same level as the ``jupyter`` directory). 
#
# 1. **gene_as.tsv**: Alternative splicing events per gene, adjusted p-value <= 0.05, fold change >= 1.5.
# 2. **all_gene_as.tsv**: all alternative splicing events
# 3. **gene_dge.tsv**: Differential gene expression, adjusted p-value <= 0.05, fold change >= 1.5
# 4. **genesWithCommonAS.tsv**: genes (as geneSymbol, the number of splicing events, and the number of tissues the event occurs in)
# 5. **Total_AS_by_chr.tsv**: Total alternative splicing events per chromosome
# 6. **Total_AS_by_geneSymbol.tsv**: Count the number of tissues in which specific genes show significant alternative splicing
# 7. **DGE_by_geneSymbol.tsv***: Most highly expressed genes by tissue
# 5. **Total_AS_by_tissue.tsv**: Count the number of significant splicing events per tissue
# 6. **Total_AS_by_splicingtype.tsv**: Count number of significant splicing events for each of the 5' alternative splicing categories
# 7. **SplicingIndex_chr.tsv**: Splicing index by chr (number of sigificant AS events per 1000 exons)

suppressMessages({
    options(warn = -1) 
    library(dplyr)
    library(Biobase)
    library(tibble)
    library(R.utils)
    library(rtracklayer)
})

# ## 1. Download all the rMATS results
#
# Each of the alternative splicing output files are downloaded here:

# ### 1.1 get released rMATS GTF annotations
#
# For each splicing type, the junctions are defined, so we have 5 specific annotated splicing specific junction ID annotation files:
#
# 1. **fromGTF.A3SS.txt**: annotations for the alternative 3' splice site junctions
# 2. **fromGTF.A5SS.txt**: annotations for the alternative 5' splice site junctions
# 3. **fromGTF.MXE.txt**: annotations for the mutually exclusive exon junctions
# 4. **fromGTF.RI.txt**: annotations for the retained introns junctions
# 5. **fromGTF.SE.txt**: annotations for the skipped exon junctions
#

# ### 2  Refined results
# We define **refined results** as (FC >= 1.5 and pVal <= 0.05) for the sex\*as_event coefficient result for the linear model

# ### 2.1 getTissueReduction

tissue_reduction_filename <- "../assets/tissues.tsv"
tissue_reduction <- read.table(tissue_reduction_filename, header=TRUE, sep="\t",
                               skipNul=FALSE, stringsAsFactors = FALSE)
colnames(tissue_reduction)  <- c("SMTSD","female","male","include","display_name")
tissue_reduction <- tissue_reduction[tissue_reduction$display_name != "n/a",]
tissue_reduction$display_name <- factor(tissue_reduction$display_name)
levels(tissue_reduction$display_name)
message("We extracted ", length(levels(tissue_reduction$display_name))," different tissues with at least 50 samples in both M & f")

# ### 2.2 Read in refined results and annotations

significant_results_dir = "../data/"
pattern = "model_B_sex_as_events_refined.csv"
files <- list.files(path = significant_results_dir, pattern = pattern)
as_types <- c("a3ss", "a5ss", "mxe", "ri", "se")
message("We extracted ", length(files), " model_B_sex_as_events_refined.csv files")

a3ss_annot <- read.table(file = "../data/fromGTF.A3SS.txt", sep = "\t", quote = "\"", header = T, stringsAsFactors = F)
a5ss_annot <- read.table(file = "../data/fromGTF.A5SS.txt", sep = "\t", quote = "\"", header = T, stringsAsFactors = F)
mxe_annot <- read.table(file = "../data/fromGTF.MXE.txt", sep = "\t", quote = "\"", header = T, stringsAsFactors = F)
ri_annot <- read.table(file = "../data/fromGTF.RI.txt", sep = "\t", quote = "\"", header = T, stringsAsFactors = F)
se_annot <- read.table(file = "../data/fromGTF.SE.txt", sep = "\t", quote = "\"", header = T, stringsAsFactors = F)

head(se_annot)

# ### 2.3 create_as_structure 
#
# This function doees an aggregation of the alternative splicing events - good for all events and the significantly expressed events.

create_as_structure <- function ( results_dir, files, all_or_das, pattern, tissue_reduction) {
    gene_as = data.frame()
    counts <- rep(NA, length(files))
    message("\nnumber of files:", paste(length(files)), collapse = "")
    for (i in 1:length(files)) {
       lines  <- read.table(file=paste0(results_dir, files[i]), 
                                     header = TRUE, sep = ",", quote = "\"'", skipNul = FALSE)
       if (dim(lines)[1] > 0) {
           event     <- as.vector(as.character(rownames(lines)))
           tissue1   <- gsub(pattern,"", files[i], fixed = TRUE)
           counts[i] <- dim(lines)[1]
           event_idx <- substring(event, regexpr("[0-9]+$", event))
           res       <- data.frame()
           if (grepl("^a3ss_", files[i])) {
               # remove the first 5 letters of the string 
               tissue2 <- substring(tissue1,6)
               idx <- match(event_idx, a3ss_annot$ID)
               res <- data.frame(GeneJunction <- event,
                              ASE          <- "A3SS", 
                              ASE_IDX      <- idx,
                              Tissue       <- tissue2,
                              counts       <- counts[i],
                              Display      <- tissue_reduction[tissue_reduction$SMTSD == tissue2, "display_name"],
                              GeneSymbol   <- a3ss_annot$geneSymbol[idx],
                              GeneID       <- a3ss_annot$GeneID[idx],
                              chr          <- a3ss_annot$chr[idx],
                              logFC        <- lines$logFC,
                              AveExpr      <- lines$AveExpr,
                              t            <- lines$t,
                              PValue       <- lines$P.Value,
                              AdjPVal      <- lines$adj.P.Val,
                              B            <- lines$B)
               colnames(res) <- c("GeneJunction","ASE","ASE_IDX","Tissue","counts","Display",
                                  "GeneSymbol","GeneID","chr","logFC","AveExpr","t","PValue","AdjPVal","B")
               gene_as <- rbind(gene_as,res)
            
           } else if (grepl("^a5ss_", files[i])) {
               # remove the first 5 letters of the string 
               tissue2 <- substring(tissue1,6)
               idx <- match(event_idx, a5ss_annot$ID)
               res <- data.frame(GeneJunction <- event,
                              ASE          <- "A5SS", 
                              ASE_IDX      <- idx,
                              Tissue       <- tissue2,
                              counts       <- counts[i],
                              Display      <- tissue_reduction[tissue_reduction$SMTSD == tissue2, "display_name"],
                              GeneSymbol   <- a5ss_annot$geneSymbol[idx],
                              GeneID       <- a5ss_annot$GeneID[idx],
                              chr          <- a5ss_annot$chr[idx],
                              logFC        <- lines$logFC,
                              AveExpr      <- lines$AveExpr,
                              t            <- lines$t,
                              PValue       <- lines$P.Value,
                              AdjPVal      <- lines$adj.P.Val,
                              B            <- lines$B)
               colnames(res) <- c("GeneJunction","ASE","ASE_IDX","Tissue","counts","Display",
                               "GeneSymbol","GeneID","chr","logFC","AveExpr","t","PValue","AdjPVal","B")
               gene_as <- rbind(gene_as,res)
           } else if (grepl("^mxe_", files[i])) {
               # remove the first 4 letters of the string 
               tissue2 <- substring(tissue1,5)
               idx <- match(event_idx, a3ss_annot$ID)
               res <- data.frame(GeneJunction <- event,
                              ASE          <- "MXE", 
                              ASE_IDX      <- idx,
                              Tissue       <- tissue2,
                              counts       <- counts[i],
                              Display      <- tissue_reduction[tissue_reduction$SMTSD == tissue2, "display_name"],
                              GeneSymbol   <- mxe_annot$geneSymbol[idx],
                              GeneID       <- mxe_annot$GeneID[idx],
                              chr          <- mxe_annot$chr[idx],
                              logFC        <- lines$logFC,
                              AveExpr      <- lines$AveExpr,
                              t            <- lines$t,
                              PValue       <- lines$P.Value,
                              AdjPVal      <- lines$adj.P.Val,
                              B            <- lines$B)
               colnames(res) <- c("GeneJunction","ASE","ASE_IDX","Tissue","counts","Display",
                                  "GeneSymbol","GeneID","chr","logFC","AveExpr","t","PValue","AdjPVal","B")
               gene_as <- rbind(gene_as,res)
           } else if (grepl("^se_", files[i])) {
               # remove the first 3 letters of the string 
               tissue2 <- substring(tissue1,4)
               idx <- match(event_idx, se_annot$ID)
               res <- data.frame(GeneJunction <- event,
                              ASE          <- "SE", 
                              ASE_IDX      <- idx,
                              Tissue       <- tissue2,
                              counts       <- counts[i],
                              Display      <- tissue_reduction[tissue_reduction$SMTSD == tissue2, "display_name"],
                              GeneSymbol   <- se_annot$geneSymbol[idx],
                              GeneID       <- se_annot$GeneID[idx],
                              chr          <- se_annot$chr[idx],
                              logFC        <- lines$logFC,
                              AveExpr      <- lines$AveExpr,
                              t            <- lines$t,
                              PValue       <- lines$P.Value,
                              AdjPVal      <- lines$adj.P.Val,
                              B            <- lines$B)
               colnames(res) <- c("GeneJunction","ASE","ASE_IDX","Tissue","counts","Display",
                                  "GeneSymbol","GeneID","chr","logFC","AveExpr","t","PValue","AdjPVal","B")
               gene_as <- rbind(gene_as,res)
           } else if (grepl("^ri_", files[i])){
               # remove the first 3 letters of the string 
               tissue2 <- substring(tissue1,4)
               idx <- match(event_idx, ri_annot$ID)
               res <- data.frame(GeneJunction <- event,
                              ASE          <- "RI", 
                              ASE_IDX      <- idx,
                              Tissue       <- tissue2,
                              counts       <- counts[i],
                              Display      <- tissue_reduction[tissue_reduction$SMTSD == tissue2, "display_name"],
                              GeneSymbol   <- ri_annot$geneSymbol[idx],
                              GeneID       <- ri_annot$GeneID[idx],
                              chr          <- ri_annot$chr[idx],
                              logFC        <- lines$logFC,
                              AveExpr      <- lines$AveExpr,
                              t            <- lines$t,
                              PValue       <- lines$P.Value,
                              AdjPVal      <- lines$adj.P.Val,
                              B            <- lines$B)
               colnames(res) <- c("GeneJunction","ASE","ASE_IDX","Tissue","counts","Display",
                                  "GeneSymbol","GeneID","chr","logFC","AveExpr","t","PValue","AdjPVal","B")
               gene_as <- rbind(gene_as,res)
           }
        
       } #if has sig. events
    
   } #for all files
   colnames(gene_as) <- c("GeneJunction","ASE","ASE_IDX","Tissue","counts","Display","GeneSymbol",
                           "GeneID","chr","logFC","AveExpr","t","PValue","AdjPVal","B")
   n_unique_genes <- length(summary(as.factor(gene_as$GeneSymbol),maxsum=50000))
   message("For the run for ", all_or_das, " run")
   message("We extracted a total of ",nrow(gene_as)," alternative splicing events (gene_as)")
   message("This includes ", n_unique_genes, " total genes")
   return (gene_as)
}

# ### 2.4 create_dge_structure 
#
# This function does an aggregation of the differential gene expression events - good for all events and the significantly expressed events.

create_dge_structure <- function ( results_dir, files, all_or_dge, pattern, map_pattern, tissue_reduction) {
   gene_dge = data.frame()
   counts <- rep(NA, length(files))
   for (i in 1:length(files)) {
      lines  <- read.table(file=paste0(results_dir, files[i]), 
                                     header = TRUE, sep = ",", quote = "\"'", skipNul = FALSE)
      if (dim(lines)[1] > 0) {
         tissue1    <- gsub(pattern,"", files[i], fixed = TRUE)
         map_lines  <- read.table(file=paste0(paste0(results_dir, tissue1),map_pattern),
                                     header = TRUE, sep = ",", quote = "\"'", skipNul = FALSE)
         counts[i]  <- dim(lines)[1]    
         tissue1    <- gsub(pattern,"", files[i], fixed = TRUE)
         map_lines  <- read.table(file=paste0(paste0(results_dir, tissue1),map_pattern),
                                     header = TRUE, sep = ",", quote = "\"'", skipNul = FALSE)
         ensg_ver   <- as.vector(as.character(rownames(lines)))
         ensg_no_ver<- as.vector(as.character(map_lines$ensg_names))
         ensg_genes <- as.vector(as.character(map_lines$ensg_genes))
         counts[i]  <- dim(lines)[1]  
         res <- data.frame(Tissue       <- tissue1,
                           ENSG_ver     <- ensg_ver,
                           ENSG_no_ver  <- ensg_no_ver,
                           GeneSymbol   <- ensg_genes,
                           counts       <- counts[i],
                           Display      <- tissue_reduction[tissue_reduction$SMTSD == tissue1, "display_name"],
                           logFC        <- lines$logFC,
                           AveExpr      <- lines$AveExpr,
                           t            <- lines$t,
                           PValue       <- lines$P.Value,
                           AdjPVal      <- lines$adj.P.Val,
                           B            <- lines$B)
         colnames(res) <- c("Tissue","ENSG_ver","ENSG_no_ver","GeneSymbol","counts","Display",
                            "logFC","AveExpr","t","PValue","AdjPVal","B")
         gene_dge <- rbind(gene_dge, res)
       } #if has sig. events
    } #for all files
    colnames(gene_dge) <- c("Tissue","ENSG_ver","ENSG_no_ver","GeneSymbol","counts","Display",
                        "logFC","AveExpr","t","PValue","AdjPVal","B")
    n_unique_genes <- length(summary(as.factor(gene_dge$GeneSymbol),maxsum=50000))
    message("For the run for ", all_or_dge, "run")
    message("We extracted a total of ",nrow(gene_dge)," gene events (gene_dge)")
    message("This includes ", n_unique_genes, " total genes")
    return(gene_dge)
}

# ### 2.5 Read in the alternative splicing results
#
# We will create an aggregation of  all the results and all the significant results
#

# +
results_dir         <- "../data/"
significant_pattern <- "_AS_model_B_sex_as_events_refined.csv"
significant_files   <- list.files(path = results_dir, pattern = significant_pattern)
all_pattern         <- "_AS_model_B_sex_as_events.csv"
all_files           <- list.files(path = results_dir, pattern = all_pattern)
as_types            <- c("a3ss", "a5ss", "mxe", "ri", "se")
message("Length of all_files: ", length(all_files))
message("Length of significant_files: ", length(significant_files))

gene_as     <- create_as_structure (results_dir      <- results_dir, 
                                    files            <- significant_files,
				    all_or_das       <- "differentially significant alternative splicing",
                                    pattern          <- significant_pattern, 
                                    tissue_reduction <- tissue_reduction)
all_gene_as <- create_as_structure (results_dir      <- results_dir, 
                                    files            <- all_files, 
				    all_or_das       <- "all alternatively spliced",
                                    pattern          <- all_pattern, 
                                    tissue_reduction <- tissue_reduction)
head(gene_as,2)
gene_as$Tissue <- factor(gene_as$Tissue)
write.table(gene_as, "../data/gene_as.tsv", quote=FALSE, sep="\t")
write.table(all_gene_as, "../data/all_gene_as.tsv", quote=FALSE, sep="\t")
# -

# ### 2.6 Create a genes-id file capturing the unique gene-junction locations in a single file
# rMATS 3.2.5 unique junction ids by splicing event tied together with gene names and these identifiers useful for downstream analyses and investigations.

results_dir         <- "../data/"
significant_pattern <- "^se_*AS_model_B_sex_as_events_refined.csv"
files   <- list.files(path =results_dir, pattern = glob2rx(significant_pattern))
message("The first file from ^se_*AS_model_B_sex_as_events_refined.csv: ", files[1])
pattern="_AS_model_B_sex_as_events_refined.csv"
geneids <- data.frame()
for (i in 1:length(files)) {
    lines  <- read.table(file=paste0(results_dir, files[i]), 
                                     header = TRUE, sep = ",", quote = "\"'", skipNul = FALSE)
    
    if (dim(lines)[1] > 0) {
           event     <- as.vector(as.character(rownames(lines)))
           tissue1   <- gsub(pattern,"", files[i], fixed = TRUE)
           event_idx <- substring(event, regexpr("[0-9]+$", event))
           res       <- data.frame()
           tissue2 <- substring(tissue1,4)
           idx <- match(event_idx, se_annot$ID)
           res <- data.frame(geneIDs      <- event,
                             ID           <- event_idx,
                             GeneSymbol   <- se_annot$geneSymbol[idx],
                             GeneID       <- se_annot$GeneID[idx],
                             chr          <- se_annot$chr[idx])
           outfilename <- paste0(paste0("../data/se_",tissue2),"_geneids.tsv")
           write.table(res, outfilename, quote=FALSE, sep="\t")
           
     }
}
message("Done writing ", length(files), " files.")

# ### 2.7 Read in the differential gene expression results
#
# Here we create an aggregation of al the significant results differential gene expression events.

# +
results_dir             <- "../data/"
significant_dge_pattern <- "_DGE_refined.csv"
significant_dge_files   <- list.files(path = results_dir, pattern = significant_dge_pattern)
map_pattern             <- "_DGE_ensg_map.csv"
length(significant_files)

gene_dge     <- create_dge_structure (results_dir      <- results_dir, 
                                      files            <- significant_dge_files, 
 				      all_or_dge       <- "differential gene expression",
                                      pattern          <- significant_dge_pattern, 
                                      map_pattern      <- map_pattern,
                                      tissue_reduction <- tissue_reduction)

head(gene_dge,2)
gene_dge$Tissue <- factor(gene_dge$Tissue)
write.table(gene_dge,     "../data/gene_dge.tsv",     quote=FALSE, sep="\t")

# Note the all_gene_dge are in the assets directory and do not change 
# write.table(all_gene_dge, "../assets/all_gene_dge.tsv", quote=FALSE, sep="\t")
# -

# Load in the gencode.v30.annotation.gtf file for additional annotation

# +
# add chr information for summary data later, use the annotation we used for rMATS
message("downloading gencode v30 annotation\n")
system("wget -O ../data/gencode.v30.annotation.gtf.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz")
message("Done!\n")
message("Unzipping compressed file gencode.v30.annotation.gtf.gz..")
system("gunzip ../data/gencode.v30.annotation.gtf.gz", intern = TRUE)
message("Done! gencode.v30.annotation.gtf can be found in ../data/")

gencode <- import("../data/gencode.v30.annotation.gtf")
gtf.df <- as.data.frame (gencode)
chr_genes <- unique(gtf.df[,c("seqnames","gene_name","gene_id")])
colnames(chr_genes) <- c("chr","GeneSymbol", "ENSG")
head(chr_genes)
# -

for (i in 1:dim(chr_genes)[1]) {
    chr_genes$ENSG[i] <- as.character(strsplit(chr_genes$ENSG[i],'\\.\\w+$'))
 #    chr_genes$ENSG[i] <- gsub('\\.*$', '', chr_genes$ENSG[i])
}
head(chr_genes)

# +
i = 1
chr <- rep("NA",dim(gene_dge)[1])
gene_dge$chr <- chr
for (i in 1:dim(gene_dge)[1]) {
    match  <- as.character(chr_genes$ENSG) %in% as.character((gene_dge$ENSG_no_ver[i]))
    if (sum(match==TRUE)== 1) {
        chr[i] <- as.character(chr_genes[match,]$chr)
        gene_dge$chr[i] <- chr[i]
    } else if (sum(match==TRUE)>1) {
        all <- as.vector(as.character(chr_genes[match,]$chr))
        gene_dge$chr[i] <- as.character(all[1])
    }
}
head(gene_dge)

write.table(gene_dge, "../data/gene_dge.tsv", quote=FALSE, sep="\t")
# -

# ### 2.8 Summary gene_as and gene_dge regarding events by tissue

# +
XY <- gene_as %>% group_by(Tissue) %>% tally()
XY <- XY[order(XY$n),decreasing=TRUE]
head(XY)
message("Minimum splicing events per tissue ", min(XY$n), " maximum splicing events per tissue ", max(XY$n))
#XY[XY$n>100,]
sum(XY$n<100)

XY <- gene_dge %>% group_by(Tissue) %>% tally()
XY <- XY[order(XY$n),decreasing=TRUE]
head(XY)
message("Minimum gene expression events per tissue ", min(XY$n), " maximum gene expression events per tissue ", max(XY$n))
# table(gtf.df[,c("gene_type")])
# -

# ### 3 Data Structures for Figures

# ### 3.1 gene_as.tsv
#
# This file contains (description)
# Here is a typical line
# <pre>
# A data.frame: 6 × 15
# GeneJunction	ASE	ASE_IDX	Tissue	counts	Display	GeneSymbol	GeneID	chr	logFC	AveExpr	t	PValue	AdjPVal	B
# <fct>	<fct>	<int>	<fct>	<int>	<fct>	<fct>	<fct>	<fct>	<dbl>	<dbl>	<dbl>	<dbl>	<dbl>	<dbl>
# 1	XIST-2253	A3SS	2253	adipose_subcutaneous	4	Adipose (sc)	XIST	ENSG00000229807.11	chrX	-4.4086049	3.196317	-36.488970	4.635568e-154	3.893877e-150	310.016049
# 2	XIST-2252	A3SS	2252	adipose_subcutaneous	4	Adipose (sc)	XIST	ENSG00000229807.11	chrX	-2.4147126	3.647690	-21.921057	1.444102e-78	6.065229e-75	160.028167
# 3	GREB1L-4933	A3SS	4933	adipose_subcutaneous	4	Adipose (sc)	GREB1L	ENSG00000141449.14	chr18	1.2793173	2.115005	7.123138	3.052112e-12	8.545914e-09	16.692429
# 4	RHCG-1776	A3SS	1776	adipose_subcutaneous	4	Adipose (sc)	RHCG	ENSG00000140519.14	chr15	-0.6930009	1.636472	-3.922124	9.797866e-05	3.919146e-02	1.142232
# 5	XIST-2253	A3SS	2253	adipose_visceral_omentum	12	Adipose (v)	XIST	ENSG00000229807.11	chrX	-4.4403352	3.113532	-33.950800	2.654474e-123	2.209585e-119	241.826117
# 6	XIST-2252	A3SS	2252	adipose_visceral_omentum	12	Adipose (v)	XIST	ENSG00000229807.11	chrX	-2.4506832	3.650617	-18.890779	2.817671e-58	1.172715e-54	114.731682
# </pre>
# There are 2887 significant events in the file.

glimpse(gene_as)
gene_as$Tissue <- factor(gene_as$Tissue)
length(levels(gene_as$Tissue))
table(is.na(gene_as$Display))
table(gene_as$Display)
colnames(gene_as)
head(gene_as)
tissue_reduction$display_name <- factor(tissue_reduction$display_name)

# +
x_as_events <- gene_as[gene_as$chr=="chrX",]
message("There were ",nrow(gene_as)," total significant alternative splicing events (gene_as)")
message("There were ",nrow(x_as_events)," total significant alternative splicing events on the X chromosome (gene_as)")
message("i.e., ", (100*nrow(x_as_events)/nrow(gene_as)), "% of all significant AS events were on the X chromosome")

numberOfUniqueTissues <- length(summary(as.factor(gene_as$Display),maxsum=500))
numberOfASEmechanisms <- length(summary(as.factor(gene_as$ASE),maxsum=500))

message("gene_as now has ",numberOfUniqueTissues, " tissues and ", numberOfASEmechanisms, " ASE categories")
message("ASE:")
summary(as.factor(gene_as$ASE),maxsum=500)
# -

# ### 3.2 gene_dge.tsv
#
# This file contains (description)
# Here is a typical line
# <pre>
# Tissue  ENSG_ver        ENSG_no_ver     GeneSymbol      counts  Display logFC   AveExpr t       PValue  AdjPVal B
# 1       adipose_subcutaneous    ENSG00000176728.7       ENSG00000176728 TTTY14  765     Adipose (sc)    -7.98216577151896     -0.928812923511535       -139.823010017733       0       0       1107.42326360464
# 2       adipose_subcutaneous    ENSG00000231535.5       ENSG00000231535 LINC00278       765     Adipose (sc)    -6.09542040758638      -2.77656379347601       -126.913818678612       0       0       1050.36559888639
# 3       adipose_subcutaneous    ENSG00000129824.15      ENSG00000129824 RPS4Y1  765     Adipose (sc)    -9.6641901864726      4.63528767282141 -125.827094717734       0       0       1041.87660796556
# 4       adipose_subcutaneous    ENSG00000067646.11      ENSG00000067646 ZFY     765     Adipose (sc)    -9.50458982938477     0.672755457406984        -125.037143030325       0       0       1033.61131617113
# 5       adipose_subcutaneous    ENSG00000229807.10      ENSG00000229807 XIST    765     Adipose (sc)    9.89280986473167      1.23756039627052 121.69689757218 0       0       1030.17757492281
# 6       adipose_subcutaneous    ENSG00000229236.1       ENSG00000229236 TTTY10  765     Adipose (sc)    -6.20901295440725     -2.74524363170072        -122.540297482165       0       0       1029.41424065532
# 7       adipose_subcutaneous    ENSG00000233864.7       ENSG00000233864 TTTY15  765     Adipose (sc)    -8.19361688496523     -0.741097276206495       -122.47199454746        0       0       1027.58175958647
# 8       adipose_subcutaneous    ENSG00000260197.1       ENSG00000260197 AC010889.1      765     Adipose (sc)    -8.52835806068555      -0.686009457030557      -119.790486729538       0       0       1015.85291786821
# 9       adipose_subcutaneous    ENSG00000183878.15      ENSG00000183878 UTY     765     Adipose (sc)    -9.52139275438866     1.60375445153084 -110.599936868261       0       0       953.069478992754
# </pre>
# There are 7417 significant events in the file.

# +
x_dge_events <- gene_dge[gene_dge$chr=="chrX",]
message("There were ",nrow(gene_dge)," total significant differential gene expression events (gene_dge)")
message("There were ",nrow(x_dge_events)," total significant differential gene expression events on the X chromosome (gene_as)")
message("i.e., ", (100*nrow(x_dge_events)/nrow(gene_dge)), "% of all significant DGE events were on the X chromosome")

numberOfUniqueTissues <- length(summary(as.factor(gene_dge$Display),maxsum=500))

message("gene_dge now has ",numberOfUniqueTissues, " tissues")
# -

# ### 3.3 Count events by chromosome
#
# Count the number of significant alternative splicing events per chromosome and save to the file **Total_AS_by_chr.tsv**.

# ### 3.3.1 by alternative splicing events

total_as_by_chr <- gene_as          %>% 
                   group_by(chr)    %>% 
                   count(chr)       %>% 
                   arrange(desc(n)) %>% 
                   as.data.frame()
total_as_by_chr$chr <- factor(total_as_by_chr$chr, levels = total_as_by_chr$chr)
length(total_as_by_chr$chr)
total_as_by_chr
glimpse(total_as_by_chr)
write.table(total_as_by_chr, file= "../data/Total_AS_by_chr.tsv", sep="\t", quote = FALSE, row.names=F)

# ### 3.3.2 by gene expression

total_dge_by_chr <- gene_dge          %>% 
                   group_by(chr)    %>% 
                   count(chr)       %>% 
                   arrange(desc(n)) %>% 
                   as.data.frame()
total_dge_by_chr$chr <- factor(total_dge_by_chr$chr, levels = total_dge_by_chr$chr)
length(total_dge_by_chr$chr)
total_dge_by_chr
glimpse(total_dge_by_chr)
write.table(total_dge_by_chr, file= "../data/Total_DGE_by_chr.tsv", sep="\t", quote = FALSE, row.names=F)

# ### 3.4 Count events by genes 

# ### 3.4.1 by alternative splicing

total_as_by_geneSymbol <- gene_as %>% 
                          group_by(GeneSymbol) %>% 
                          count(GeneSymbol)    %>% 
                          arrange(desc(n))     %>% 
                          as.data.frame()
total_as_by_geneSymbol$GeneSymbol <- factor(total_as_by_geneSymbol$GeneSymbol, 
                                            levels = total_as_by_geneSymbol$GeneSymbol)
length(total_as_by_geneSymbol$GeneSymbol)
head(total_as_by_geneSymbol,10)
write.table(total_as_by_geneSymbol, file = "../data/Total_AS_by_geneSymbol.tsv", sep = "\t", quote=FALSE, row.names = F)
# ### 3.5 Count most frequent splicing by tissue

# ### 3.5.1 by alternative splicing

total_as_by_tissue <- gene_as %>% 
                      group_by(Display) %>% 
                      count(Display)    %>% 
                      arrange(desc(n))  %>% 
                      as.data.frame()
total_as_by_tissue$Display <- factor(total_as_by_tissue$Display, 
                                     levels = total_as_by_tissue$Display)
head(total_as_by_tissue,10)
length(total_as_by_tissue$Display)
write.table(total_as_by_tissue, file = "../data/Total_AS_by_tissue.tsv", sep = "\t", row.names = F)

# ### 3.5.2 by gene expression

#glimpse(gene_dge)
gene_dge$GeneSymbol <- factor(gene_dge$GeneSymbol)
total_dge_by_tissue <- gene_dge %>% 
                          select(c(GeneSymbol, Display, logFC)) %>%
                          group_by(Display) %>%
                          arrange(desc(logFC)) %>%
                          tally() %>%
                          arrange(desc(n)) %>%
                          as.data.frame()
head(total_dge_by_tissue,10)
length(total_dge_by_tissue$Display)
write.table(total_dge_by_tissue, file = "../data/Total_DGE_by_tissue.tsv", sep = "\t", quote=FALSE, row.names = F)

# ###  3.6 Significant Count by splicing type 
# We define **significant** to be FC > 1.5 and pVal < 0.05
#
# Our starting values were the significant events, all meeting the criteria FC > 1.5 and pVal < 0.05
#

total_as_by_splicingtype <- gene_as %>% 
                            group_by(ASE)    %>% 
                            count(ASE)       %>% 
                            arrange(desc(n)) %>%
                            as.data.frame()
total_as_by_splicingtype$ASE <- factor(total_as_by_splicingtype$ASE, levels = total_as_by_splicingtype$ASE)
total_as_by_splicingtype
write.table(total_as_by_splicingtype, file= "../data/Total_AS_by_splicingtype.tsv")

# ###  3.7 Significant Count by splicing type (significant == FC > 1.5 and pVal < 0.05)

# +
A3SS_keep <- as.character(gene_as$ASE) %in% "A3SS"
table(A3SS_keep)
A3SS.gene_as <- data.frame(gene_as[A3SS_keep == TRUE,])

A5SS_keep <- as.character(gene_as$ASE) %in% "A5SS"
table(A5SS_keep)
A5SS.gene_as <- data.frame(gene_as[A5SS_keep == TRUE,])

MXE_keep  <- as.character(gene_as$ASE) %in% "MXE"
table(MXE_keep)
MXE.gene_as <- data.frame(gene_as[MXE_keep == TRUE,])

SE_keep   <- as.character(gene_as$ASE) %in% "SE"
table(SE_keep)
SE.gene_as <- data.frame(gene_as[SE_keep == TRUE,])

RI_keep   <- as.character(gene_as$ASE) %in% "RI"
table(RI_keep)
RI.gene_as <- data.frame(gene_as[RI_keep == TRUE,])

dim(A3SS.gene_as)
dim(A5SS.gene_as)
dim(MXE.gene_as)
dim(SE.gene_as)
dim(RI.gene_as)

# -
# ### 3.8 Siginficant spliced by Gene for each splicing factor

# +
A3SS.res <- A3SS.gene_as %>% group_by(GeneSymbol) %>% count(GeneSymbol) %>% arrange(desc(n)) %>% as.data.frame()
A3SS.res$GeneSymbol <- factor(A3SS.res$GeneSymbol, levels = A3SS.res$GeneSymbol)
message("Significant spliced genes for A3SS\n",
        paste(length(A3SS.res$GeneSymbol)), collapse=" ")
#head(A3SS.res)

A5SS.res <- A5SS.gene_as %>% group_by(GeneSymbol) %>% count(GeneSymbol) %>% arrange(desc(n)) %>% as.data.frame()
A5SS.res$GeneSymbol <- factor(A5SS.res$GeneSymbol, levels = A5SS.res$GeneSymbol)
message("Significant spliced genes for A5SS\n",
        paste(length(A5SS.res$GeneSymbol)), collapse=" ")
#head(A5SS.res)

MXE.res <- MXE.gene_as %>% group_by(GeneSymbol) %>% count(GeneSymbol) %>% arrange(desc(n)) %>% as.data.frame()
MXE.res$GeneSymbol <- factor(MXE.res$GeneSymbol, levels = MXE.res$GeneSymbol)
message("Significant spliced genes for MXE\n",
        paste(length(MXE.res$GeneSymbol)), collapse=" ")
#head(MXE.res)

RI.res <- RI.gene_as %>% group_by(GeneSymbol) %>% count(GeneSymbol) %>% arrange(desc(n)) %>% as.data.frame()
RI.res$GeneSymbol <- factor(RI.res$GeneSymbol, levels = RI.res$GeneSymbol)
message("Significant spliced genes for RI\n",
        paste(length(RI.res$GeneSymbol)), collapse=" ")
#head(RI.res)

SE.res <- SE.gene_as %>% group_by(GeneSymbol) %>% count(GeneSymbol) %>% arrange(desc(n)) %>% as.data.frame()
SE.res$GeneSymbol <- factor(SE.res$GeneSymbol, levels = SE.res$GeneSymbol)
message("Significant spliced genes for SE\n",
        paste(length(SE.res$GeneSymbol)), collapse=" ")
#head(SE.res)
# -

# ### 3.9 Count most frequent spliced genes

# +
genesMostFrequentlySpliced <- gene_as %>% 
                              group_by(GeneSymbol) %>% 
                              count(GeneSymbol)    %>% 
                              arrange(desc(n))     %>% 
                              as.data.frame()
genesMostFrequentlySpliced$GeneSymbol <- factor(genesMostFrequentlySpliced$GeneSymbol, 
                                                levels = genesMostFrequentlySpliced$GeneSymbol)
length(genesMostFrequentlySpliced$GeneSymbol)

#Add number of tissues
nTissues <- rep(NA, length(genesMostFrequentlySpliced))
for (i in 1:nrow(genesMostFrequentlySpliced)) {
    df_gene <- gene_as %>% 
               filter(GeneSymbol == genesMostFrequentlySpliced$GeneSymbol[i])
    nTissues[i] <- length(unique(df_gene$Tissue))
}
genesMostFrequentlySpliced$Tissues <- nTissues
head(genesMostFrequentlySpliced)
write.table(genesMostFrequentlySpliced, file = "../data/genesWithCommonAS.tsv", sep = "\t", quote = F, row.names = F)
# -

# ### 3.10 Count most frequent spliced chromosomes
# To get an indication of which chromosome has the most frequent slicing event (regardless of type)
# We create an index based upon the number of exons per chromosome.
#
# get the annotation file, at this writing, gencode.v30.annotation.gtf
# The information as to the number of exons within the chromosome may be found there

# +
exons <- gencode[ gencode$type == "exon", ]
exons <- as.data.frame(exons)

#Obtain chromosomes we have splicing information for (recall we did not use chr Y in our analysis)
all_chr <- as.character(unique(gene_as$chr))
chr_counts <- rep(0, length(all_chr))


for (i in 1:length(all_chr)) {
  chr_counts[i] <- nrow(exons[exons$seqnames == all_chr[i], ])
}

exon_counts <- data.frame(chr = all_chr, counts = chr_counts)

# Count most frequent spliced chromosomes
res <- gene_as %>% group_by(chr) %>% count(chr) %>% arrange(desc(n)) %>% as.data.frame()
res$chr <- factor(res$chr, levels = res$chr)

idx <- match(res$chr, exon_counts$chr)

res$ExonCounts <- exon_counts$counts[idx]

res$Index <- (res$n / res$ExonCounts) * 1000

res_sorted <- res %>% arrange(desc(Index))
res_sorted$chr <- factor(res_sorted$chr, levels = res_sorted$chr)
glimpse(res_sorted)
write.table(res_sorted, file = "../data/SplicingIndex_chr.tsv", sep = "\t", quote = F, row.names = F)
# -

# ### 3.11 Overlap between Differential Gene Expression and Differential Alternative Splicing
#
# First gather the data

total_AS_Genes <- read.table(file="../data/Total_AS_by_geneSymbol.tsv", header=TRUE, sep="\t",
                               skipNul=FALSE, stringsAsFactors = FALSE)
sigAsGenes <- sort(total_AS_Genes$GeneSymbol)
dge <- read.table("../data/gene_dge.tsv", sep = "\t", header = FALSE, row.names=1, skip = 1)
dge_genes <- sort(dge$V5)
head(dge_genes)
all_genes_data <- read.table("../assets/all_gene_dge.tsv")
names(all_genes_data) <- c("GeneSymbol", "ensg")
all_genes <- sort(all_genes_data$GeneSymbol)

head(all_genes)

# ### 3.12 We then do a hypergeometric/Fisher test to look for overrepresentation
# The universe consists of all genes with at least one read (all_genes_data).
# So we have
#
# |  	|  DGE+| DGE-|
# |-	|-	|-	|
# | DAS+|  a|  b|
# | DAS-|  c| d|

message("Number of sigAsGenes ", length(sigAsGenes))
notSigAs <- setdiff(all_genes,sigAsGenes)
message("Number of genes that are NOT sigAs ", length(notSigAs))
message("Number of DGE genes", length(dge_genes))
notDGE <- setdiff(all_genes,dge_genes)
message("Number of genes that are NOT DGE ", length(notDGE))
a <- intersect(sigAsGenes, dge_genes)
b <- intersect(sigAsGenes, notDGE)
c <- intersect(notSigAs, dge_genes)
d <- intersect(notSigAs, notDGE)
message("a: ", length(a), "; b: ",  length(b), "; c: ",  length(c), "; d: ",  length(d))

m <- matrix(c(length(a),length(b),length(c),length(d)), nrow=2,byrow = TRUE)
fisher.test(m)

# ### Appendix - Metadata
#
# For replicability and reproducibility purposes, we also print the following metadata:
#
# 1. Checksums of **'artefacts'**, files generated during the analysis and stored in the folder directory **`data`**
# 2. List of environment metadata, dependencies, versions of libraries using `utils::sessionInfo()` and [`devtools::session_info()`](https://devtools.r-lib.org/reference/session_info.html)

# ### Appendix 1. Checksums with the sha256 algorithm

# +
rm (notebookid)
notebookid   = "countGenesAndEvents"
notebookid

message("Generating sha256 checksums of the file `../data/gene_as.tsv` directory .. ")
system(paste0("cd ../data && find . -name gene_as.tsv -exec sha256sum {} \\;  >  ../metadata/", notebookid, "_sha256sums.txt"), intern = TRUE)
message("Done!\n")

message("Generating sha256 checksums of the file `../data/all_gene_as.tsv` directory .. ")
system(paste0("cd ../data && find . -name gene_dge.tsv -exec sha256sum {} \\;  >  ../metadata/", notebookid, "_sha256sums.txt"), intern = TRUE)
message("Done!\n")

message("Generating sha256 checksums of the file `../data/gene_dge.tsv` directory .. ")
system(paste0("cd ../data && find . -name gene_as.tsv -exec sha256sum {} \\;  >  ../metadata/", notebookid, "_sha256sums.txt"), intern = TRUE)
message("Done!\n")

message("Generating sha256 checksums of the file `../data/Total_AS_by_chr.tsv` directory .. ")
system(paste0("cd ../data && find . -name Total_AS_by_chr.tsv -exec sha256sum {} \\;  >  ../metadata/", notebookid, "_sha256sums.txt"), intern = TRUE)
message("Done!\n")

message("Generating sha256 checksums of the file `../data/Total_AS_by_geneSymbol.tsv` directory .. ")
system(paste0("cd ../data && find . -name Total_AS_by_geneSymbol.tsv -exec sha256sum {} \\;  >  ../metadata/", notebookid, "_sha256sums.txt"), intern = TRUE)
message("Done!\n")

message("Generating sha256 checksums of the file `../data/Total_AS_by_tissue.tsv` directory .. ")
system(paste0("cd ../data && find . -name Total_AS_by_tissue.tsv -exec sha256sum {} \\;  >  ../metadata/", notebookid, "_sha256sums.txt"), intern = TRUE)
message("Done!\n")

message("Generating sha256 checksums of the file `../data/Total_AS_by_splicingtype.tsv` directory .. ")
system(paste0("cd ../data && find . -name Total_AS_by_splicingtype.tsv -exec sha256sum {} \\;  >  ../metadata/", notebookid, "_sha256sums.txt"), intern = TRUE)
message("Done!\n")

message("Generating sha256 checksums of the file `../data/genesWithCommonAS.tsv` directory .. ")
system(paste0("cd ../data && find . -name genesWithCommonAS.tsv -exec sha256sum {} \\;  >  ../metadata/", notebookid, "_sha256sums.txt"), intern = TRUE)
message("Done!\n")

message("Generating sha256 checksums of the file `../data/SplicingIndex_chr.tsv` directory .. ")
system(paste0("cd ../data && find . -name SplicingIndex_chr.tsv -exec sha256sum {} \\;  >  ../metadata/", notebookid, "_sha256sums.txt"), intern = TRUE)
message("Done!\n")

# -

# ### Appendix 2. Libraries metadata

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


