# Getting set up fast on the Robinson-Lab Jupyter Lab notebook.

## 1 Start your jupyter notebook

After getting your notebook started, do the following in a terminal window
Conda install is your friend.   Saying yes at the end means you can fire and forget.

### 1.1 These conda install commands below get you going to start. 

conda install -c bioconda bioconductor-edger -y
conda install -c conda-forge r-readr -y
conda install -c bioconda bioconductor-multtest -y
conda install -c conda-forge r-r.utils -y 
conda install -c conda-forge r-snakecase -y
conda install -c conda-forge r-statmod -y
conda install -c bioconda bioconductor-stan -y
conda install -c bioconda bioconductor-rtracklayer -y
conda install -c conda-forge r-reshape -y

## 1.2 What do you do if there are more missing packages

Refer to the documentation to go further with Christina’s gitbook.

https://lifebit.gitbook.io/cloudos/web-interface/jupyter-lab#how-can-i-install-the-dependencies-i-need-in-the-notebook-session

# 2 Getting the files

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








