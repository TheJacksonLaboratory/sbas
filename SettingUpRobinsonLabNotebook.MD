# Getting set up fast on the Robinson-Lab Jupyter Lab notebook.

## 1 Start your jupyter notebook

After getting your notebook started, do the following in a terminal window
Conda install is your friend.   Saying yes at the end means you can fire and forget.

### 1.1 These conda install commands below get you going to start. 

```bash
conda install -c bioconda bioconductor-edger -y
conda install -c conda-forge r-readr -y
conda install -c bioconda bioconductor-multtest -y
conda install -c conda-forge r-r.utils -y 
conda install -c conda-forge r-snakecase -y
conda install -c conda-forge r-statmod -y
conda install -c bioconda bioconductor-stan -y
conda install -c bioconda bioconductor-rtracklayer -y
conda install -c conda-forge r-reshape -y
```

## 1.2 What do you do if there are more missing packages

Refer to the documentation to go further with Christina’s gitbook.

https://lifebit.gitbook.io/cloudos/web-interface/jupyter-lab#how-can-i-install-the-dependencies-i-need-in-the-notebook-session

# 2 Cloning the repository

You should begin by git cloning the repository

## 2.1 git clone the Robinson Lab sbas repository


# 3 Getting the files

You can get the data from the jupyternotebooks if you use piggyback or we have the in the following buckets, in the cloud environment.
To access the results from differential gene expression analysis and differential splicing junction analysis.


The results from rmats-nf NextFlow pipeline analysis of the GTEx V8 release are found in github release repositories but for convenience of our analysis now are in the google bucket named robinson-bucket in the folder named notebooks, organised in .tar.gz files


There are 5 files.

| filename | Description  |
| ------------ | ------------- |
| data.tar.gz | Contains results from the rmats-nf pipeline and the results from the execution of the first differential gene expression pipeline as captured in sbas/jupyter/differentialGeneExpressionAnalysis.ipynb |
| assets.tar.gz | This file contains the tissues.tsv file that was used for generating display labels (see the code for how to do this in the current sbas/jupyter/figure1.ipynb notebook - soon to be named sexBiasedASEventsPerTissue.ipynb). This reduces our tissue list to 39 based upon number of samples with adequate data.  Selection is indicated by the column `include` where 1 indicates keep and 0 indicates exclude |
| DGE_splicing_data.tar.gz | Results from the execution of the differential splicing junction analysis as executed for each of the reduced 39 tissue as outlined in assets.tar.gz.   The file details are explained in the next section |
| DGE_splicing_pdf.tar.gz | These are the MDS_plots -- useful for investigating the ability of voom to clean up the segregation of results -- a summary figure will be included in the publication of some subset of this information. |

## 3.1 Contents of data.tar.gz

Also available through github using wget or piggyback within R - location noted in box


| filename | Description  |
| ------------ | ------------- |
| DGE_gene_csv.tar.gz | Results from jupyter/differentialGeneExpression.ipynb |
| fromGTF.tar.gz | Annotation files from rMATS 3.2.5 on the gencode v.30 GTF file |
| gtex.corrected.rds | Corrected Yarn generated expressionSet object on the gene expression data from GTEx v8 with the phenotype data available - |
| rmats_final.a3ss.jc.ijc.txt.gz | Inclusion junction counts for the a3ss |
| rmats_final.a3ss.jc.inclen.txt.gz | rMATS percent spliced in for a3ss |
| rmats_final.a3ss.jc.sjc.txt.gz | Skipped junction counts for a3ss |
| rmats_final.a3ss.jc.skiplen.txt.gz | Length of junction for the alternative a3ss |
| rmats_final.a5ss.jc.ijc.txt.gz | Inclusion junction counts for the a5ss  |
| rmats_final.a5ss.jc.inclen.txt.gz | Inclusion junction length for a5ss |
| rmats_final.a5ss.jc.inc.txt.gz | rMATS percent spliced in for a5ss |
| rmats_final.a5ss.jc.sjc.txt.gz | Skipped junction counts for a5ss |
| rmats_final.a5ss.jc.skiplen.txt.gz | Skipped junction length for a5ss |
| rmats_final.mxe.jc.ijc.txt.gz | Inclusion junction counts for mxe |
| rmats_final.mxe.jc.inclen.txt.gz | length of junction for the mxe |
| rmats_final.mxe.jc.inc.txt.gz | rMATS percent spliced in for mxe |
| rmats_final.mxe.jc.sjc.txt.gz | Skipped junction counts mxe |
| rmats_final.mxe.jc.skiplen.txt.gz | Skipped junction length for mxe |
| rmats_final.ri.jc.ijc.txt.gz | Inclusion junction counts for ri |
| rmats_final.ri.jc.inclen.txt.gz | Inclusion junction length for ri |
| rmats_final.ri.jc.inc.txt.gz | rMATS percent spliced in for ri |
| rmats_final.ri.jc.sjc.txt.gz |  Skipped junction counts for ri |
| rmats_final.ri.jc.skiplen.txt.gz | Skipped junction length |
| rmats_final.se.jc.ijc.txt.gz | Inclusion junction counts for se |
| rmats_final.se.jc.inclen.txt.gz | Inclusion junction length for se |
| rmats_final.se.jc.inc.txt.gz | rMATS percent spliced in for se |
| rmats_final.se.jc.sjc.txt.gz | Skipped junction counts for se |
| rmats_final.se.jc.skiplen.txt.gz | Skipped junction length |
| SraRunTable.txt.gz | Annotation file connecting SRR with SAMPID |
| srr_pdata.csv | Updated phenodata matching SAMPID to SRR (1 to many) |
| Tissue_list.csv | Listing of tissues for this study in assets put as part of this release for convenience |
